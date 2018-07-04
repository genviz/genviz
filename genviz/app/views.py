import json
import pickle
import random
import textwrap
from django.core import serializers
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.generic import TemplateView, View
from Bio import Entrez
from Bio import SeqIO
from .models import *
import xmltodict
import hgvs.parser
from django.contrib.auth import login, authenticate
from django.shortcuts import render, redirect
from genviz.forms import *
from django.shortcuts import render, get_object_or_404
from .utils.variations import *

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

            handle_summary = Entrez.esummary(id=','.join(record['IdList']), db='nucleotide', rettype='gb', retmode='text')
            
            res = list(Entrez.parse(handle_summary))
            return self.render_to_response(context={
                'results': res,
                'gene': gene,
                'organism': organism
            })

        return self.render_to_response(context)

class GeneSearch(TemplateView):
    template_name = 'search.html'

class GeneDetails(TemplateView):
    template_name = 'details.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        seq_id = request.GET.get('id', None)

        if seq_id is not None:
            Entrez.email = "email@example.com"
            
            handle = Entrez.efetch(id=seq_id, db='nucleotide', rettype='gb', retmode='text')
            res = SeqIO.read(handle, format='gb')

            acc_id = res.name
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

            # Fetch variations from Clinvar
            clinvar_variations = fetch_clinvar_variations(res.id)
            dbsnp_variations = fetch_dbsnp_variations(res.id)
            external_variations = []
            for v in clinvar_variations + dbsnp_variations:
                # TODO: Support all operations
                if v.operation == 'unknown':
                    continue
                if v.coordinate_type == 'c' and 'CDS' in features:
                    cds_start = features['CDS'][0]['location'][0]
                    v.start += cds_start
                    v.end += cds_start
                if v.start >= 0 and v.end <= gene_length:
                    external_variations.append(v)

            # Get variations
            user_variations = list(Variation.objects.filter(acc_id=acc_id).all())
            variations = external_variations + user_variations
            return self.render_to_response(context={
                'entry': res,
                'entry_dict': res.__dict__,
                'features_json': json.dumps(features),
                'features': features,
                'sequence': str(res.seq),
                'gene_length': gene_length,
                'seq_id': seq_id,
                'variation_sources': set(map(lambda a: a.source, variations)),
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
    variations = Variation.objects.all()
    return render(request,'profile.html',{'patients':patients,'variations':variations})

def var(request):
    variations = Variation.objects.all()
    return render(request,'var.html',{'variations':variations})

def predict(request, pathology_id):
    pathology = Pathology.objects.get(pk=pathology_id)
    model = pickle.load(open(pathology.prediction_model, 'rb'))
    prediction = model.predict([[random.randint(0, 1) for i in range(model.n_features_)]])[0]
    return JsonResponse({
        'prediction': str(prediction),
        'precision': str(0.5723)
    })

def patient_detail(request, pk):
    patient = get_object_or_404(Patient, pk=pk)
    variations = Variation.objects.filter(patient=patient)
    pathologies = Pathology.objects
    #anot = Variation.objects.filter(patient=patient)
    return render(request, 'patient_detail.html', {
        'patient': patient,
        'variations': variations,
        'pathologies': pathologies
    })