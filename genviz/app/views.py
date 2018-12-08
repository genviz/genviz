import hgvs.parser
import json
import os
import pickle , csv
import random
import re
import textwrap
import xmltodict
import requests, sys
from Bio import Entrez
from Bio import SeqIO
from django.core import serializers
from django.core.urlresolvers import reverse_lazy, reverse
from django.http import (
    HttpResponse, 
    JsonResponse, 
    HttpResponseRedirect,
    HttpResponseNotFound
)
from django.forms.models import modelform_factory
from django.views.generic import RedirectView, TemplateView, ListView, View
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
import urllib.request


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

class TcgaSearchResults(TemplateView):
    template_name = 'tcga_results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        gene = request.GET.get('gene', None)

        if gene:
            contents = requests.get(
                "http://exac.hms.harvard.edu/rest/awesome?query=" + gene,
                headers={
                    "Content-Type": "application/json"
                }
            )
            
            if not contents.ok:
                return HttpResponseNotFound('<h1>Gene Not Found</h1>')
            
            contents = contents.json()
            transcripts_in_gene = contents['transcripts_in_gene'] 
            coverage_stats = contents['coverage_stats']
            gene_info = contents['gene']
            transcript = contents['transcript']
            variations = contents['variants_in_transcript']
            return self.render_to_response(context={
               "transcripts_in_gene" : transcripts_in_gene,
               "gene_info" : gene_info,
               "transcript": transcript,
               "variations": variations,
               "gene": gene
            })

        return self.render_to_response(context)

class EnsemblSearchResults(TemplateView):
    template_name = 'ensembl_results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        gene = request.GET.get('gene', None)
        server = "http://exac.hms.harvard.edu/rest/awesome?query="

        if gene:

            r = requests.get(
                server + gene, 
                headers={
                    "Content-Type": "application/json"
                }
            )

            if not r.ok:
                return HttpResponseNotFound('<h1>Gene Not Found</h1>')

            decoded = r.json()
            return self.render_to_response({
                "gene": gene,
                "ensembl_response": decoded["variants_in_transcript"]
            })

        return self.render_to_response(context)


class EnsemblVariantSearchResults(TemplateView):
    template_name = 'ensembl_var_results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        variations = request.GET.getlist('id')
        print(variations)
        server = "https://rest.ensembl.org"

        if variations:
            response = []
            samples = []
            ext = "/variation/human/"
            params = "?genotypes=1;content-type=application/json"

            for var in variations:
                r = requests.get(server + ext + var + params)
                decoded = r.json()
                tmp = SampleName(list(map(lambda d: d['sample'], decoded['genotypes'])))
                samples = Sample.objects.filter(individual_id__in=tmp).values('individual_id', 'population', 'gender')
                decoded['genotypes'] = samples
                response.append(decoded)

            
            return self.render_to_response({
                "response": response
            })

        return self.render_to_response(context)

def SampleName(samples):
    names = []
    for sample in samples:
        tmp = sample.split(":")
        tmp[-1] = tmp[-1].replace(" ","_").replace(".","")
        names.append(tmp[-1])
    return names

class GeneSearch(TemplateView):
    template_name = 'search.html'

class SnpSearch(TemplateView):
    template_name = 'search_snp.html'


class GeneDetailsByName(RedirectView):

    def get(self, request, *args, **kwargs):

        gene = kwargs['gene']
        print (gene)
        organism = request.GET.get('organism', None) or 'Homo sapiens'

        if gene and organism:
            Entrez.email = "email@example.com"

            term = 'RefSeqGene[keyword] AND "{}"[gene] AND "{}"[orgn]'.format(
                gene, organism)
            handle = Entrez.esearch(term=term, db='nucleotide', idtype='acc')
            record = Entrez.read(handle)
            ids = record['IdList']
            handle_summary = Entrez.esummary(id=','.join(
                ids), db='nucleotide', rettype='gb', retmode='text')

            try:
                print("Entro en try")
                res = list(Entrez.parse(handle_summary))
                print("Logro parse en try", res[0])
            except RuntimeError:
                res = []

            gene_id = res[0]['Id']

            return HttpResponseRedirect(reverse('details') + '?id=' + gene_id)


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

class UploadView(TemplateView):
    template_name = 'upload.html'

    def post(self, request, *args, **kwargs):
        name = request.POST.get('name', None)
        file = request.FILES['file']

        model = WFile.objects.create(name=name, file=file, user=request.user)
 
        return HttpResponseRedirect(request.POST.get('next', '/'))

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
    _fields = (
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
    _widgets = {
        'birthday': forms.DateInput(attrs={'class':'datepicker'}),
        'sample_date': forms.DateInput(attrs={'class':'datepicker'}),
    }
    form_class = modelform_factory(Patient, fields=_fields, widgets=_widgets)

    def form_valid(self, form):
        """
        If the form is valid, save the associated model.
        """
        self.object = form.save()
        self.request.user.patients.add(self.object)
        return super(PatientsNew, self).form_valid(form)

class PatientsUpdate(UpdateView):
    model = Patient
    template_name = 'patients/patient_form.html'
    success_url = reverse_lazy('patients_list')
    _fields = (
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
    _widgets = {
        'birthday': forms.DateInput(attrs={'class':'datepicker'}),
        'sample_date': forms.DateInput(attrs={'class':'datepicker'}),
    }
    form_class = modelform_factory(Patient, fields=_fields, widgets=_widgets)

class PatientsDelete(DeleteView):
    model = Patient
    template_name = 'patients/patient_confirm_delete.html'
    success_url = reverse_lazy('patients_list')

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


class SampleDetails(TemplateView):
    template_name = 'sample_details.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        individual_id = kwargs['string']
        try:
            sample = Sample.objects.get(individual_id=individual_id)
            context['sample'] = sample
        except:
            pass

        return self.render_to_response(context)
