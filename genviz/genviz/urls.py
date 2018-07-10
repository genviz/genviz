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
    url(r'^details/', GeneDetails.as_view(), name='details'),
    url(r'^variations/', VariationsView.as_view()),
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^patient/new/$', views.patient_new, name='patient_new'),
    url(r'^profile/$', views.profile ,name='profile'),
    url(r'^var/$', views.var ,name='var'),
    url(r'^predict/(?P<pathology_id>[\d]+)$', views.predict ,name='predict'),
    url(r'^patient/(?P<patient_id>[\w]+)/$', views.patient_detail, name='patient_detail'),
]



