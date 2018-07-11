import textwrap
import re
from Bio import Entrez
from ..models import *
import xmltodict
import hgvs.parser
import collections

Entrez.email = 'example@example.com'

def snp_to_hgvs(snp, acc_id=None):
    """
    Fetch SNP and return a list of the variations related to it as Variation objects

    @param    snp: str|[str]  SNP ID(s)
    @param acc_id: str  If provided, get variations only on specified Accession ID
    
    @return [Variation]
    """
    handle = Entrez.efetch(id=snp, db='snp', retmode='xml')
    dbsnp_xml = handle.read()
    dbsnp_dict = xmltodict.parse(dbsnp_xml)
    hgvsparser = hgvs.parser.Parser()
    variations = []
    variation = dbsnp_dict['ExchangeSet']['Rs']
    for hgvs_var in variation['hgvs']:
        if not acc_id or acc_id in hgvs_var:
            try:
                variations.append(hgvsparser.parse_hgvs_variant(hgvs_var))
            except Exception as e:
                print('(dbSNP) Couldn\'t parse variation {}'.format(hgvs_var))
                print(repr(e))
    return variations

def feature_name(feature):
    """
    Given a feature, return a friendly display name

    @param feature: dict  Feature from parsed sequence

    @return str  Feature name
    """
    if 'number' in feature['qualifiers']:
        return '{feature_type} {number}'.format(feature_type=feature['type'], number=feature['qualifiers']['number'][0])
    elif 'standard_name' in feature['qualifiers']:
        return '{feature_type} ({name})'.format(feature_type=feature['type'], name=feature['qualifiers']['standard_name'][0])
    return feature['type']

def fetch_clinvar_variations(acc_id):
    """
    Fetch ClinVar variations associated to an Accession

    @param acc_id: str  Accession ID of a sequence

    @return [Variation]
    """
    handle = Entrez.esearch(term='%s[Nucleotide/Protein Accession]' % acc_id, db='clinvar', retmax='1000')
    res = Entrez.read(handle, validate=False)
    var_ids = res['IdList']
    handle = Entrez.efetch(id=var_ids, db='clinvar', rettype='variation')
    clinvar_xml = handle.read()
    clinvar_dict = xmltodict.parse(clinvar_xml)
    hgvsparser = hgvs.parser.Parser()
    variations = []
    if 'ClinVarResult-Set' in clinvar_dict:
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
                        print('(Clinvar) Couldn\'t parse variation {}'.format(v['#text']))
                        print(repr(e))
    return variations

def fetch_dbsnp_variations(acc_id, start=1, end=1e12):
    """
    Fetch dbSNP variations associated to an Accession

    @param acc_id: str  Accession ID of a sequence
    @param  start: int  Start position to look for variations
    @param    end: int  End position to look for variations

    @return [Variation]
    """

    # If it's a chromosome
    if acc_id.startswith('NC'):
        chromosome = re.sub(r'^NC_0+(\d{1,2})(?:\.\d+)', r'\1', acc_id)
        term = '%s[Chromosome] AND (%s[CHRPOS]:%s[CHRPOS])' % (chromosome, start, end)
    else:
        term = '"%s"[Accession] AND (pathogenic[Clinical_Significance] OR other[Clinical_Significance] OR "likely benign"[Clinical_Significance] OR "likely pathogenic"[Clinical_Significance] OR benign[Clinical_Significance])' % acc_id
    print("Query on dbSNP: ", term)
    handle = Entrez.esearch(term=term, db='snp', retmax='1000')
    res = Entrez.read(handle, validate=False)
    var_ids = res['IdList']
    print('{} SNP fetched'.format(len(var_ids)))
    return fetch_snp(var_ids, acc_id) if var_ids else []

def fetch_snp(snp, acc_id=None):
    """
    Given a SNP, fetch it and return the list of variations associated to it
    with extra information such as URL and comment

    @param    snp: str  SNP ID
    @param acc_id: str  If provided, get variations only on specified Accession ID

    @return [Variation]
    """
    # TODO: Refactor, reuse snp_to_hgvs
    try:
        handle = Entrez.efetch(id=snp, db='snp', retmode='xml')
        dbsnp_xml = handle.read()
        dbsnp_dict = xmltodict.parse(dbsnp_xml)
        hgvsparser = hgvs.parser.Parser()
        variations = []
        # Handle cases with one and with multiple results
        res = dbsnp_dict['ExchangeSet']['Rs']
        rs_list = res if type(res) == list else [res]
        for variation in rs_list:
            for hgvs_var in variation['hgvs']:
                if (not acc_id) or acc_id in hgvs_var:
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
    except Exception as e:
        return []
