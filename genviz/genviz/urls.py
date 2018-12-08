"""

genviz URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
    
"""

from django.conf.urls import url, include
from django.contrib import admin
from app.views import *
from app import views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^accounts/', include('django.contrib.auth.urls')),
    url(r'^search/$', GeneSearch.as_view(), name='search'),
    url(r'^search/snp/$', SnpSearch.as_view()),
    url(r'^/?$', Home.as_view()),
    url(r'^results/$', GeneSearchResults.as_view(), name='gene_results'),
    url(r'^results/snp/$', SnpSearchResults.as_view(), name='snp_results'),
    url(r'^results/tcga/$', TcgaSearchResults.as_view(), name='tcga_results'),
    url(r'^results/1000-genomes/$', EnsemblSearchResults.as_view(), name='ensembl_results'),
    url(r'^results/1000-genomes-variations/$', EnsemblVariantSearchResults.as_view(), name='ensembl_variant_results'),
    url(r'^goto/details/(?P<gene>\w+)$', GeneDetailsByName.as_view(), name='details_by_name'),
    url(r'^details/', GeneDetails.as_view(), name='details'),
    url(r'^variations/', VariationsView.as_view()),
    url(r'^upload/', UploadView.as_view(), name='upload'),
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^predict/(?P<pathology_id>[\d]+)$', views.predict ,name='predict'),
    url(r'^patients/new/$', PatientsNew.as_view(), name='patients_new'),
    url(r'^patients/$', PatientsList.as_view(), name='patients_list'),
    url(r'^patients/(?P<pk>\d+)$', PatientsDetail.as_view(), name='patients_detail'),
    url(r'^patients/edit/(?P<pk>\d+)$', PatientsUpdate.as_view(), name='patients_edit'),
    url(r'^patients/delete/(?P<pk>\d+)$', PatientsDelete.as_view(), name='patients_delete'),
    url(r'^sample/details/(?P<string>\w+)$', SampleDetails.as_view(), name='sample_details'),
]



