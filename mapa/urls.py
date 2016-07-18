from django.conf.urls import url, include #dodato include
from .import views
from django.views.generic import TemplateView #Dodato


urlpatterns = [
	
    #url(r'^/а$',  views.Mapa.as_view(), name='index'),
    url(r'^$', 'mapa.views.mapa_pocetna'),
    #url(r'^stanice$', TemplateView.as_view(template_name='stanice.html')),
    #url(r'^stanice.json$', 'mapa.views.dajStanice'),
    					

]