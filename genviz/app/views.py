import json
from django.core import serializers
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.generic import TemplateView, View
from Bio import Entrez
from Bio import SeqIO
from .models import *

# TODO: Create view to show details of a gene
def fetch_gene_details(ids):
    handle_genes = Entrez.efetch(id=','.join(ids), db='nucleotide', rettype='gb', retmode='text')
    records_genes = SeqIO.parse(handle_genes, 'gb')
    res = {r.id: str(r.seq) for r in records_genes}

class GeneSearchResults(TemplateView):
    template_name = 'results.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        gene = request.GET.get('gene', None)
        organism = request.GET.get('organism', None) or 'Homo sapiens'

        if gene and organism:
            Entrez.email = "email@example.com"

            term = '"{}"[gene] AND "{}"[orgn]'.format(gene, organism)        
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

            # Add coordinates to plot gene features
            features = {}
            for f in res.features:
                f_type = f.type
                f = f.__dict__
                f['location'] = (f['location'].start, f['location'].end)
                features[f_type] = features.setdefault(f_type, []) + [f]

            for f_type in features:
                features[f_type].sort(key=lambda f: f['location'][1])
            
            gene_length = features['source'][0]['location'][1]
            sequence = []   
            prev_end = 0
            if 'CDS' in features:
                for feature in features['CDS']:
                    cds_start, cds_end = feature['location']
                    translation = feature['qualifiers']['translation'][0]
                    sequence.append({
                        'type': 'non-CDS',
                        'location': [int(prev_end), int(cds_start)],
                        'sequence': str(res.seq[prev_end:cds_start])
                    })
                    sequence.append({
                        'type': 'CDS',
                        'location': list(map(int, feature['location'])),
                        'sequence': str(res.seq[cds_start:cds_end]),
                        'translation': feature['qualifiers']['translation'][0],
                        'triplets': []
                    })
                    for i in range(len(translation)):
                        sequence[-1]['triplets'].append({
                            'sequence': str(res.seq[cds_start+i*3:cds_start+(i+1)*3]),
                            'translation': feature['qualifiers']['translation'][0][i]
                        })
                    # If there's a leftover base at the end
                    if len(sequence[-1]['triplets'])*3+1 == cds_end-cds_start+1:
                        sequence[-1]['triplets'].append({
                            'sequence': str(res.seq[cds_end]),
                            'translation': '-'    
                        })

                    prev_end = cds_end

            if prev_end != gene_length:
                sequence.append({
                    'type': 'non-CDS',
                    'location': [int(prev_end), int(gene_length)],
                    'sequence': str(res.seq)
                })

            # Get annotations
            annotations = Annotation.objects.filter(seq_id=seq_id).all()
            return self.render_to_response(context={
                'entry': res,
                'entry_dict': res.__dict__,
                'features_json': json.dumps(features),
                'features': features,
                'sequence': sequence,
                'sequence_json': json.dumps(sequence),
                'gene_length': gene_length,
                'seq_id': seq_id,
                'annotations_json': json.dumps(list(map(lambda x: x["fields"], json.loads(serializers.serialize('json', annotations)))))
            })

class AnnotationsView(View):
    def post(self, request, *args, **kwargs):
        annotations_json = json.loads(request.POST.get('annotations', '[]'))
        seq_id = request.POST.get('seq_id', None)
        for annotation_json in annotations_json:
            annotation = Annotation(author=request.user, seq_id=seq_id, **annotation_json)
            annotation.save()
        return HttpResponseRedirect(request.POST.get('next', '/'))
