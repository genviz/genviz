import hgvs.parser
import json
import os
import pickle
import random
import re
import textwrap
import xmltodict
from Bio import Entrez
from Bio import SeqIO
from django.core import serializers
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.generic import TemplateView, ListView, View
from django.views.generic.detail import DetailView
from django.views.generic.edit import (
    CreateView,
    UpdateView,
    DeleteView
)
from django.contrib.auth import login, authenticate
from django.shortcuts import render, redirect, get_object_or_404
from genviz.forms import *
from genviz.settings import BASE_DIR
from .models import *
from .utils.variations import *
from .utils.prediction import *


class Home(TemplateView):
    template_name = 'home.html'

class GeneSearchResults(TemplateView):
    template_name = 'results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        gene = request.GET.get('gene', None)
        organism = request.GET.get('organism', None) or 'Homo sapiens'

        if gene and organism:
            Entrez.email = "email@example.com"

            term = 'RefSeqGene[keyword] AND "{}"[gene] AND "{}"[orgn]'.format(gene, organism)        
            handle = Entrez.esearch(term=term, db='nucleotide', idtype='acc')
            record = Entrez.read(handle)
            ids = record['IdList']
            handle_summary = Entrez.esummary(id=','.join(ids), db='nucleotide', rettype='gb', retmode='text')
            
            other_term = '"transcript variant"[All fields] AND "{}"[gene] AND "{}"[orgn]'.format(gene, organism)        
            other_handle = Entrez.esearch(term=other_term, db='nucleotide', idtype='acc')
            other_record = Entrez.read(other_handle)
            other_ids = list(set(other_record['IdList']) - set(ids))
            other_handle_summary = Entrez.esummary(id=','.join(other_ids), db='nucleotide', rettype='gb', retmode='text')
    
            try:
                res = list(Entrez.parse(handle_summary))
            except RuntimeError:
                res = []

            try:
                other_res = sorted(
                    list(filter(
                        lambda r: 'transcript variant' in r['Title'].lower(),
                        Entrez.parse(other_handle_summary
                    ))),
                    key=lambda r: r['Title']
                )
            except RuntimeError:
                other_res = []

            return self.render_to_response(context={
                'results': res,
                'other_results': other_res,
                'gene': gene,
                'organism': organism
            })

        return self.render_to_response(context)


class SnpSearchResults(TemplateView):
    template_name = 'snp_results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        snp = request.GET.get('snp', None)

        if snp:
            Entrez.email = "email@example.com"
            variations = fetch_snp(snp)
            acc_ids = set(map(lambda v: v.acc_id, variations))
            # Get name of sequences
            handle = Entrez.esummary(id=','.join(acc_ids), db='nucleotide', rettype='gb', retmode='text')
            res = list(Entrez.parse(handle))
            acc_ids_titles = { x['AccessionVersion']: x['Title'] for x in res }
            return self.render_to_response(context={
                'variations': variations,
                'snp': snp,
                'acc_ids_titles': acc_ids_titles
            })

        return self.render_to_response(context)


class GeneSearch(TemplateView):
    template_name = 'search.html'

class SnpSearch(TemplateView):
    template_name = 'search_snp.html'

class GeneDetails(TemplateView):
    template_name = 'details.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        seq_id = request.GET.get('id', None)
        start = request.GET.get('start', None)
        end = request.GET.get('end', None)

        if seq_id is not None:
            Entrez.email = "email@example.com"
            params = {
                'id': seq_id,
                'db': 'nucleotide',
                'rettype': 'gbwithparts',
                'retmode': 'text'
            }
            if start and end:
                params.update({
                    'seq_start': start,
                    'seq_stop': end    
                })
            handle = Entrez.efetch(**params)
            res = SeqIO.read(handle, format='gb')
            acc_id = res.id
            # Add coordinates to plot gene features
            features = {}
            for f in res.features:
                f_type = f.type
                f = f.__dict__
                f['location'] = (f['location'].start, f['location'].end)
                f['display_name'] = feature_name(f)
                features[f_type] = features.setdefault(f_type, []) + [f]

            for f_type in features:
                features[f_type].sort(key=lambda f: f['location'][1])

            gene_length = features['source'][0]['location'][1] - features['source'][0]['location'][0] + 1
            start = max(int(start), 1) if start else 1
            end = int(end) if end else gene_length

            # Fetch variations from Clinvar
            clinvar_variations = fetch_clinvar_variations(res.id)
            print('Loading dbSNP variations...')
            dbsnp_variations = fetch_dbsnp_variations(res.id, start, end)
            print('Done dbSNP variations')
            external_variations = []
            variations_databases = set()
            for v in clinvar_variations + dbsnp_variations:
                # TODO: Support all operations
                if v.operation == 'unknown':
                    continue
                if v.coordinate_type == 'c' and 'CDS' in features:
                    cds_start = features['CDS'][0]['location'][0]
                    v.start += cds_start
                    v.end += cds_start
                if v.start >= start and v.end <= end:
                    external_variations.append(v)
                    variations_databases.add(v.source)

            # Get variations
            user_variations = list(Variation.objects.filter(acc_id=acc_id).all())
            variations = external_variations + user_variations

            print('Rendering template')
            return self.render_to_response(context={
                'entry': res,
                'start': start,
                'end': end,
                'entry_dict': res.__dict__,
                'features_json': json.dumps(features),
                'features': features,
                'sequence': str(res.seq),
                'gene_length': gene_length,
                'seq_id': seq_id,
                'variations_databases': variations_databases,
                'variations_patients': set(map(lambda a: a.source, user_variations)),
                'variations_json': json.dumps([x["fields"] for x in json.loads(serializers.serialize('json', variations))]),
                'patients': request.user.patients if request.user.is_authenticated and request.user.is_biologist else []
           })

class VariationsView(View):
    def post(self, request, *args, **kwargs):
        variations_json = json.loads(request.POST.get('variations', '[]'))
        acc_id = request.POST.get('acc_id', None)
        for variation_json in variations_json:
            patient = get_object_or_404(Patient, pk=variation_json['patient_id'])
            variation = Variation(
                author=request.user,
                source="({}) {}".format(patient.identifier, patient.full_name()),
                acc_id=acc_id,
                **variation_json
            )
            variation.save()
        return HttpResponseRedirect(request.POST.get('next', '/'))


def signup(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.save()
            #falta login
            return redirect('/accounts/login/')
    else:
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})

    
def patient_new(request):
    if request.method == 'POST':
        form = PatientForm(request.POST)
        if form.is_valid():
            patient = form.save()
            request.user.patients.add(patient)
            request.user.save()
            return redirect('/profile')
    else:
        form = PatientForm()
    return render(request, 'patient/patient_new.html', {'form': form})

def profile(request):
    patients = request.user.patients.all()
    patients_json = json.dumps([x["fields"] for x in json.loads(serializers.serialize('json', patients))])
    variations = Variation.objects.all()
    return render(
        request,
        'profile.html',
        {
            'patients': patients,
            'patients_json': patients_json,
            'variations': variations
        }
    )

def var(request):
    variations = Variation.objects.all()
    return render(request,'var.html',{'variations':variations})

def predict(request, pathology_id):
    pathology = Pathology.objects.get(pk=pathology_id)
    patient_id = request.GET['patient_id']
    patient = Patient.objects.get(pk=patient_id)
    path = os.path.join(BASE_DIR, pathology.prediction_model)
    model = pickle.load(open(path, 'rb'))
    vector = feature_from_patient(pathology, patient)
    print("Predicting with ", vector)
    prediction = model.predict([vector])[0]
    if prediction == 0:
        precision = pathology.precision_negative
    else:
        precision = pathology.precision_positive
    return JsonResponse({
        'prediction': str(prediction),
        'precision': str(precision)
    })

def patient_detail(request, patient_id):
    patient = get_object_or_404(Patient, identifier=patient_id)
    variations = Variation.objects.filter(patient=patient)
    pathologies = Pathology.objects
    #anot = Variation.objects.filter(patient=patient)
    return render(request, 'patient_detail.html', {
        'patient': patient,
        'variations': variations,
        'pathologies': pathologies
    })

class PatientsList(ListView):
    model = Patient
    template_name = 'patients/patient_list.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        patients = self.request.user.patients.all()
        patients_json = json.loads(serializers.serialize('json', patients))
        patients_json = json.dumps([dict({'pk': x['pk']}, **x["fields"]) for x in patients_json])

        context['patients_json'] = patients_json
        return context

class PatientsDetail(DetailView):
    model = Patient
    template_name = 'patients/patient_detail.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        variations = Variation.objects.filter(patient=context['patient'])
        pathologies = Pathology.objects

        context['pathologies'] = pathologies
        context['variations'] = variations
        return context


class PatientsNew(CreateView):
    model = Patient
    template_name = 'patients/patient_new.html'
    success_url = reverse_lazy('patients_list')
    fields = (
        'identifier',
        'first_name',
        'last_name',
        'sex',
        'birthday',
        'sample_date',
        'phone',
        'diagnosis',
        'email'
    )
    widgets = {
        'birthday': forms.DateInput(attrs={'class':'datepicker'}),
        'sample_date': forms.DateInput(attrs={'class':'datepicker'}),
    }

class PatientsUpdate(UpdateView):
    model = Patient
    success_url = reverse_lazy('patients_list')
    fields = (
        'identifier',
        'first_name',
        'last_name',
        'sex',
        'birthday',
        'sample_date',
        'phone',
        'diagnosis',
        'email'
    )
    widgets = {
        'birthday': forms.DateInput(attrs={'class':'datepicker'}),
        'sample_date': forms.DateInput(attrs={'class':'datepicker'}),
    }

class PatientsDelete(DeleteView):
    model = Patient
    success_url = reverse_lazy('patients_list')