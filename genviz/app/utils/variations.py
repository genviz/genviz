import textwrap
from Bio import Entrez
from ..models import *
import xmltodict
import hgvs.parser

Entrez.email = 'example@example.com'

def snp_to_hgvs(snp):
    handle = Entrez.efetch(id=snp, db='snp', retmode='xml')
    dbsnp_xml = handle.read()
    dbsnp_dict = xmltodict.parse(dbsnp_xml)
    hgvsparser = hgvs.parser.Parser()
    variations = []
    variation = dbsnp_dict['ExchangeSet']['Rs']
    for hgvs_var in variation['hgvs']:
        if acc_id in hgvs_var:
            try:
                variations.append(hgvsparser.parse_hgvs_variant(hgvs_var))
            except Exception as e:
                print('Couldn\'t parse variation {}'.format(hgvs_var))
                print(repr(e))
    return variations

def feature_name(feature):
    """
    Given a feature, return a friendly display name
    """
    if 'number' in feature['qualifiers']:
        return '{feature_type} {number}'.format(feature_type=feature['type'], number=feature['qualifiers']['number'][0])
    elif 'standard_name' in feature['qualifiers']:
        return '{feature_type} ({name})'.format(feature_type=feature['type'], name=feature['qualifiers']['standard_name'][0])
    return feature['type']

def fetch_gene_details(ids):
    handle_genes = Entrez.efetch(id=','.join(ids), db='nucleotide', rettype='gb', retmode='text')
    records_genes = SeqIO.parse(handle_genes, 'gb')
    res = { r.id: str(r.seq) for r in records_genes }

def fetch_clinvar_variations(acc_id):
    handle = Entrez.esearch(term='%s[Nucleotide/Protein Accession]' % acc_id, db='clinvar')
    res = Entrez.read(handle, validate=False)
    var_ids = res['IdList']
    handle = Entrez.efetch(id=var_ids, db='clinvar', rettype='variation')
    clinvar_xml = handle.read()
    clinvar_dict = xmltodict.parse(clinvar_xml)
    hgvsparser = hgvs.parser.Parser()

    variations = []
    for variation in clinvar_dict['ClinVarResult-Set']['VariationReport']:
        for v in variation['Allele']['HGVSlist']['HGVS']:
            if '@AccessionVersion' in v and v['@AccessionVersion'] == acc_id:
                try:
                    clinical_significance = variation['ClinicalAssertionList']['GermlineList']['Germline']['ClinicalSignificance']['Description']
                    comment = textwrap.dedent("""\
                    {hgvs}
                    Clinical significance: {significance}
                    
                    (click for more information)
                    """.format(hgvs=v['#text'], significance=clinical_significance))
                    
                    variation_obj = Variation.from_hgvs_obj(hgvsparser.parse_hgvs_variant(v['#text']), 'clinvar')
                    variation_obj.url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/{}'.format(variation['@VariationID'])
                    variation_obj.comment = comment
                    variations.append(variation_obj)
                except Exception as e:
                    print('Couldn\'t parse variation {}'.format(v['#text']))
                    print(repr(e))
    return variations

def fetch_dbsnp_variations(acc_id):
    handle = Entrez.esearch(term='%s[Nucleotide/Protein Accession]' % acc_id, db='snp')
    res = Entrez.read(handle, validate=False)
    var_ids = res['IdList']
    handle = Entrez.efetch(id=var_ids, db='snp', retmode='xml')
    dbsnp_xml = handle.read()
    dbsnp_dict = xmltodict.parse(dbsnp_xml)
    hgvsparser = hgvs.parser.Parser()
    variations = []
    for variation in dbsnp_dict['ExchangeSet']['Rs']:
        for hgvs_var in variation['hgvs']:
            if acc_id in hgvs_var:
                try:
                    #clinical_significance = variation['ClinicalAssertionList']['GermlineList']['Germline']['ClinicalSignificance']['Description']
                    comment = textwrap.dedent("""\
                    {hgvs}
                    (click for more information)
                    """.format(hgvs=hgvs_var))
                    
                    variation_obj = Variation.from_hgvs_obj(hgvsparser.parse_hgvs_variant(hgvs_var), 'dbSNP')
                    variation_obj.url = 'https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={}'.format(variation['@rsId'])
                    variation_obj.comment = comment
                    variations.append(variation_obj)
                except Exception as e:
                    print('Couldn\'t parse variation {}'.format(hgvs_var))
                    print(repr(e))
    return variations