from django.http import HttpResponse, JsonResponse
from django.views.generic import TemplateView
from Bio import Entrez
from Bio import SeqIO

# TODO: Create view to show details of a gene
def fetch_gene_details(ids):
    handle_genes = Entrez.efetch(id=','.join(ids), db='nucleotide', rettype='gb', retmode='text')
    records_genes = SeqIO.parse(handle_genes, 'gb')
    res = {r.id: str(r.seq) for r in records_genes}

class GeneSearch(TemplateView):
    template_name = 'search.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        gene = request.GET.get('gene', None)
        orgn = request.GET.get('organism', None) or 'Homo sapiens'

        if gene is not None:
            Entrez.email = "email@example.com"
            
            search_term = '"{}"[gene] AND "{}"[orgn]'.format(gene, orgn)
            handle = Entrez.esearch(term=search_term, db='nucleotide', idtype='acc')
            record = Entrez.read(handle)

            handle_summary = Entrez.esummary(id=','.join(record['IdList']), db='nucleotide', rettype='gb', retmode='text')
            
            res = list(Entrez.parse(handle_summary))
            return JsonResponse(res, safe=False)

        return self.render_to_response(context)

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

            return self.render_to_response(context={
                'entry': res,
                'entry_dict': res.__dict__
            })