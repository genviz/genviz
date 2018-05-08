from django.http import HttpResponse, JsonResponse
from django.views.generic import TemplateView
from Bio import Entrez
from Bio import SeqIO

class SearchGene(TemplateView):
    template_name = 'search.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        #import pdb; pdb.set_trace()
        term = request.GET.get('term', None)
        if term is not None:
            Entrez.email = "email@example.com"
            
            handle = Entrez.esearch(term='{}[gene]'.format(term), db='nucleotide', idtype='acc')
            record = Entrez.read(handle)
            ids = ','.join(record['IdList'][:5])
            handle_genes = Entrez.efetch(id=ids, db='nucleotide', rettype='gb', retmode='text')
            records_genes = SeqIO.parse(handle_genes, 'gb')

            res = {r.id: str(r.seq) for r in records_genes}
            return JsonResponse(res)

        return self.render_to_response(context)