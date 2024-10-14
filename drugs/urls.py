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

    url(r'^NewDrugsBrowser', views.NewDrugsBrowser.as_view(), name='NewDrugsBrowser'),
    ## These three can be the same view maybe? ##
    url(r'^drug_selection', views.DrugSectionSelection.as_view(), name='DrugSearch'),
    # url(r'^target_selection', views.DrugSectionSelection.as_view(page='targets'), name='TargetSearch'),
    # url(r'^disease_selection', views.DrugSectionSelection.as_view(page='diseases'), name='DiseaseSearch'),
    #############################################
    # url(r'^overview', views.DruggedGPCRome.as_view(), name='DruggedGPCRome'),
    ## Make a single view for the Venn diagrams, add a variable to define which page is shown ##
    url(r'^drugs_venn', views.DrugsVenn, name='DrugsVenn'),
    url(r'^targets_venn', views.TargetVenn, name='DrugsVenn'),
    ############################################################################################
    url(r'^TargetSelectionTool', views.TargetSelectionTool.as_view(), name='TargetSelectionTool')
]
