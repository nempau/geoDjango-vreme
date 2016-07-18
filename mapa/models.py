from django.db import models
#dodato kao u aplikaciji
from django.contrib.gis.db import models as gismodels


class WeatherStation(gismodels.Model):

    wmoid = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=256)

    geom = gismodels.PointField()

    objects = gismodels.GeoManager()

    def __str__(self):
        return self.name
# Create your models here.
