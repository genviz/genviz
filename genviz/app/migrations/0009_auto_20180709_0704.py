# -*- coding: utf-8 -*-
# Generated by Django 1.11.13 on 2018-07-09 07:04
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('app', '0008_auto_20180709_0649'),
    ]

    operations = [
        migrations.RenameField(
            model_name='pathology',
            old_name='model_precision_negative',
            new_name='precision_negative',
        ),
        migrations.RenameField(
            model_name='pathology',
            old_name='model_precision_positive',
            new_name='precision_positive',
        ),
    ]
