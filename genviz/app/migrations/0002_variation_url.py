# -*- coding: utf-8 -*-
# Generated by Django 1.11.13 on 2018-06-13 04:06
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('app', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='variation',
            name='url',
            field=models.CharField(max_length=255, null=True),
        ),
    ]
