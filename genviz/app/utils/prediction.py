import pandas as pd
import re
import random
from app.models import *
from datetime import datetime, timedelta
from django.db.models import Max
from .variations import *

SEX_TO_NUM = {
	'F': 0,
	'M': 1
}
NUM_TO_SEX = {v: k for k, v in SEX_TO_NUM.items()}

GENOTYPE_TO_NUM = {
	'homozygous-ref': 0,
	'heterozygous': 1,
	'homozygous': 2,
}
NUM_TO_GENOTYPE = {v: k for k, v in GENOTYPE_TO_NUM.items()}

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
	rs_ids = list(map(lambda v: re.sub(r'(rs\d+).*', r'\1', v).strip(), data[0:0]))
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

def create_patients_from_csv(csv, pathology, doctor, rows):
	"""
	@param csv: path to csv
	@param pathology: Pathology object
	@param doctor: User object
	@param rows: array with index of rows to add
	"""
	data = pd.read_csv(csv, sep=',')
	n = pathology.variations.all().aggregate(Max('header_order'))['header_order__max']

	for i in rows:
		row = data[i:i+1]
		age = int(row['Edad'])
		days_per_year = 365.24
		sample_date = datetime.now()
		patient = Patient(
			first_name = 'Test patient',
			last_name = str(i),
			identifier = 'TEST{}'.format(i),
			sex = NUM_TO_SEX[int(row['Sexo'])],
			birthday = sample_date - timedelta(days=(age*days_per_year + 100)),
			sample_date = sample_date
		)
		patient.save()
		doctor.patients.add(patient)
		# TODO: Take it out of the loop (query var_loc only once)
		# TODO: Get correct ref, alt, variation from data (derive it from genotype and rs)
		for i in range(n+1):
			var_loc = pathology.variations.filter(header_order=i)[0]
			variation = Variation(
				start=var_loc.start,
				end=var_loc.end,
				acc_id=var_loc.acc_id,
				# Column by index
				genotype=NUM_TO_GENOTYPE[int(row.iloc[:, 4+i])],
				ref=random.choice(['A', 'T', 'G', 'C']),
				alt=random.choice(['A', 'T', 'G', 'C']),
				operation='sub',
				source='TESTDATA',
				patient=patient
			)
			variation.save()
	doctor.save()


def feature_from_patient(pathology, patient):
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
	for i in range(n+1):
		if i in genotypes:
			feature.append(GENOTYPE_TO_NUM[genotypes[i]])
		else:
			feature.append(GENOTYPE_TO_NUM['homozygous-ref'])
	return feature
