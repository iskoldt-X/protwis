from django.conf.urls import url
from django.views.decorators.cache import cache_page

from drugs import views

urlpatterns = [
    url(r'^drugbrowser', cache_page(60 * 60 * 24 * 28)(views.DrugBrowser.as_view()), name='drugbrowser'),
    url(r'^drugstatistics',  views.drugstatistics, name='drugstatistics'),
    url(r'^drugmapping',  views.drugmapping, name='drugmapping'),
    url(r'^nhs/section/(?P<slug>[\w|\W]+)/$',  views.nhs_section, name='nhs_section'),
    url(r'^nhs/(?P<slug>[\w|\W]+)/$',  views.nhs_drug, name='nhs_drug'),
    url(r'^indication/(?P<code>[\w|\W]+)/$',  views.indication_detail, name='indication'),
    url(r'^DrugsIndicationsTargets', views.Drugs_Indications_Targets.as_view(), name='Drugs_Indications_Targets'),
    url(r'^Drugs', views.Drugs.as_view(), name='Drugs'),
    url(r'^Indications', views.Indications.as_view(), name='Indications'),
    url(r'^Targets', views.Targets.as_view(), name='Targets'),
    url(r'^TargetSelectionTool', views.TargetSelectionTool.as_view(), name='TargetSelectionTool')
]
