"""vreme URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Import the include() function: from django.conf.urls import url, include
    3. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf import settings
from django.conf.urls import url, include
from django.conf.urls.static import static
from django.contrib import admin
#dodejemo i tacke
from djgeojson.views import GeoJSONLayerView
from mapa.models import WeatherStation

from django.views.generic import TemplateView

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^mapa', include('mapa.urls')),
    url(r'^data.geojson$', GeoJSONLayerView.as_view(model=WeatherStation), name='data'),
    #url(r'^data/{z}/{x}/{y}.geojson$', GeoJSONLayerView.as_view(model=WeatherStation), name='data'),
    url(r'^stanice$', TemplateView.as_view(template_name='stanice.html')),
    url(r'^index2$', TemplateView.as_view(template_name='index2.html')),
    url(r'^blanko$', TemplateView.as_view(template_name='blanko.html')),
    url(r'^stanice.json$', 'mapa.views.dajStanice'),

    ]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)




#url(r'^data.geojson$', GeoJSONLayerView.as_view(model=WeatherStation), name='data'),