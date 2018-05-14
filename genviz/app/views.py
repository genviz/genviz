import json
from django.http import HttpResponse, JsonResponse
from django.views.generic import TemplateView
from Bio import Entrez
from Bio import SeqIO

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
        #import pdb; pdb.set_trace()
        gene_id = request.GET.get('id', None)

        if gene_id is not None:
            Entrez.email = "email@example.com"
            
            handle = Entrez.efetch(id=gene_id, db='nucleotide', rettype='gb', retmode='text')
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
                        'location': [prev_end, cds_start],
                        'sequence': res.seq[prev_end:cds_start] 
                    })
                    sequence.append({
                        'type': 'CDS',
                        'location': feature['location'],
                        'sequence': res.seq[cds_start:cds_end],
                        'translation': feature['qualifiers']['translation'][0],
                        'triplets': []
                    })
                    for i in range(len(translation)):
                        sequence[-1]['triplets'].append({
                            'sequence': res.seq[cds_start+i*3:cds_start+(i+1)*3],
                            'translation': feature['qualifiers']['translation'][0][i]
                        })
                    # If there's a leftover base at the end
                    if len(sequence[-1]['triplets'])*3+1 == cds_end-cds_start+1:
                        sequence[-1]['triplets'].append({
                            'sequence': res.seq[cds_end],
                            'translation': '-'    
                        })

                    prev_end = cds_end

            if prev_end != gene_length:
                sequence.append({
                    'type': 'non-CDS',
                    'location': [prev_end, gene_length],
                    'sequence': res.seq
                })

            return self.render_to_response(context={
                'entry': res,
                'entry_dict': res.__dict__,
                'features_json': json.dumps(features),
                'features': features,
                'sequence': sequence
            })