from django.http import HttpResponse
from django.shortcuts import render

from django.views.generic import TemplateView
#Pokusaj renderovanja
from mapa.models import WeatherStation
from django.contrib.gis.geos import Polygon
import json



def mapa_pocetna(request):
	#return HttpResponse('<h1> Zdravo ja sam mapa!</h1>')
	return render(request, "index.html", {})
# Create your views heea.

#Funkcija za renderovanje dela tacakatype
def dajStanice(request):
	# Pronaći minimalni obuhvatni pravougaonik

	bbox= request.GET['bbox'].split(',')
	poly=Polygon.from_bbox(bbox)

	#Povuci podatke iz baze

	stanice= WeatherStation.objects.filter(geom__within=poly)


	#Konvertuj u geoJson format
	geojson_dict= {

		"type":"FeatureCollection",
		"features":[stanice_u_geojson(stanica) for stanica in stanice],
	}
	# vrati odgovor
	return HttpResponse(json.dumps(geojson_dict),
							content_type='application/json',)

def stanice_u_geojson(stanica):
	return {
		"type":"Feature",
		"geometry": {
			#"coordinates": [5.883333333333333, 36.81666666666667],
			"type": "Point",

			"coordinates":[stanica.geom.x, stanica.geom.y ],

		},

		"properties": {
			"description": stanica.name
		},

		"id": stanica.wmoid,
	}

	'''def borjStanica(WeatherStation):
		broj= WeatherStation.objects.count()
		return broj'''

'''class Mapa(TemplateView):
    template_name = 'vreme/index.html'


# Create your views here.'''
