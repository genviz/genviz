import pandas as pd
import re
from app.models import *
from django.db.models import Max
from .variations import *

def create_pathology_from_csv(csv, name, prediction_model, precision_p, precision_n):
	"""
	@param csv:  path to csv
	@param name: pathology name
	"""
	p = Pathology(
		name = name, 
		prediction_model = prediction_model,
		precision_positive = precision_p,
		precision_negative = precision_n
	)
	p.save()
	data = pd.read_csv(csv, sep=',')
	data.drop('Subject', axis=1, inplace=True)
	data.drop('Grupo', axis=1, inplace=True)
	data.drop('HISTOPATOLÓGICO', axis=1, inplace=True)
	data.drop('Edad', axis=1, inplace=True)
	data.drop('Sexo', axis=1, inplace=True)
	# Get rs ids from headers
	rs_ids = list(map(lambda v: re.sub(r'(rs\d+).*', r'\1', v), data[0:0]))
	for i, rs_id in enumerate(rs_ids):
		variations = snp_to_hgvs(rs_id)
		for variation in variations:
			var_obj = Variation.from_hgvs_obj(variation, 'dbSNP')
			v_location = VariationLocation(
				start=var_obj.start,
				end=var_obj.end,
				acc_id=var_obj.acc_id,
				header_order=i
			)
			v_location.save()
			p.variations.add(v_location)
	p.save()

def feature_from_patient(pathology, patient):
	SEX_TO_NUM = {
		'F': 0,
		'M': 1
	}
	GENOTYPE_TO_NUM = {
		'homozygous-ref': 0,
		'heterozygous': 1,
		'homozygous': 2,
	}
	feature = [patient.age_at_sample_date(), SEX_TO_NUM[patient.sex]]
	genotypes = {}

	for v in pathology.variations.order_by('header_order').all():
		variation = Variation.objects.filter(
			patient=patient,
			start=v.start,
			end=v.end,
			acc_id=v.acc_id
		)
		if variation:
			genotypes[v.header_order] = variation[0].genotype

	n = pathology.variations.all().aggregate(Max('header_order'))['header_order__max']
	for i in range(n):
		if i in genotypes:
			feature.append(GENOTYPE_TO_NUM[genotypes[i]])
		else:
			feature.append(GENOTYPE_TO_NUM['homozygous-ref'])
	return feature
