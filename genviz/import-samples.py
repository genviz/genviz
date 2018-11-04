# open file & create csvreader
from app.models import *
import csv
# import the relevant model

with open("../samples.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        sample = Sample(
            family_id = row[0],
            individual_id = row[1],
            paternal_id = row[2],
            maternal_id = row[3],
            gender = row[4],
            population = row[6],
            relationship = row[7],
            siblings = row[8],
            second_order = row[9],
            third_order = row[10],
            comments = row[11]
        )

        sample.save()

