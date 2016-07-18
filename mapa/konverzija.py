import csv
from django.contrib.gis.geos import Point

from mapa.models import WeatherStation
from django.contrib.gis.utils import LayerMapping


def dms2dec(value):
    """
    Degres Minutes Seconds to Decimal degres
    """
    degres, minutes, seconds = value.split()
    seconds, direction = seconds[:-1], seconds[-1]
    dec = float(degres) + float(minutes)/60 + float(seconds)/3600
    if direction in ('S', 'W'):
        return -dec
    return dec

def run():
	csv_file = '/home/nemanja/Fakultet/GisProgramiranje/vreme/vreme/Pub9volA160314x.flatfile.txt'
	i=0
	reader = csv.DictReader(open(csv_file, 'rt'), delimiter="\t")
	
	for line in reader:
		lng = dms2dec(line.pop('Longitude'))
		lat = dms2dec(line.pop('Latitude'))
		wmoid = int(line.pop('StationId'))
		name = line.pop('StationName').title()

		WeatherStation (wmoid=wmoid, name=name, geom=Point(lng, lat)).save()
		if i==9500:
			break
		i=i+1
		print(i)



'''def run():
	with open('/home/nemanja/Fakultet/GisProgramiranje/vreme/mapa/vreme.txt') as csvfile:
		reader = csv.DictReader(csvfile)
		for line in reader:
			lng = dms2dec(line['Longitude'])
			lat = dms2dec(line['Latitude'])
			wmoid = int(line['StationId'])
			name = line['StationName'].title()
			WeatherStation(wmoid=wmoid, name=name, geom=Point(lng, lat)).save()'''


