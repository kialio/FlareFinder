#!/usr/bin/env python

__version__ = "0.0.1"
__author__ = "Jeremy. S. Perkins"


import numpy as np
#import asciitable as atable
from urllib import urlretrieve,urlopen
from datetime import datetime
import re
from shutil import copy2
import pyfits
from BeautifulSoup import BeautifulSoup
from astropy.io import ascii as atable
from astropy.table import Table

class FlareFinder:

	def __init__(self, verbose = False, instrument = 'LATAP', mjd = 0):

		self.verbose = verbose
		self.instrument = instrument
		self.baseURL = {"LATAP":"http://fermi.gsfc.nasa.gov/ssc/data/access/lat/2yr_catalog/ap_lcs/",
						"BAT":"http://swift.gsfc.nasa.gov/results/transients/"}


		if(mjd == 0):
			self.mjd = self.getMJD()
		else:
			self.mjd = mjd

	def getMJD(self):
	
		dt_now = datetime.now()
		dt_mjdref = datetime(1858,11,17)
	
		diff = dt_now - dt_mjdref
	
		if self.verbose:
			print "The MJD today is {}.".format(diff.days)
		
		return diff.days

	def getCurve(self,sourceName = "2FGLJ0007.7+6825c", web = True, directory = './ap_lcs', save = False):


		if self.instrument == "LATAP":

			fileName = self.baseURL['LATAP'] + "lc_" + re.sub('\+', 'p', sourceName) + ".dmp1.out"
			tmpFileName = directory + "/lc_" + re.sub('\+', 'p', sourceName) + ".dmp1.out"
			outFileName = "lc_" + re.sub('\+', 'p', sourceName) + ".dmp1.out"

		elif self.instrument == "BAT":

			fileName = self.baseURL["BAT"] + sourceName + ".lc.txt"
			tmpFileName = directory + "/" + sourceName + ".lc.txt"
			outFileName = sourceName + ".lc.txt"

		else:
			print "Unknown instrument."
			return -1

		if web:

			tmpFileName,response =  urlretrieve(fileName)
			if self.verbose:
				print fileName,response.gettype()

			if response.gettype() != 'text/plain':
				print "Didn't get a TEXT file."
				return -1 

		if(save):
			outFileName = "lc_" + re.sub('\+', 'p', sourceName) + ".dmp1.out"
			copy2(tmpFileName,outFileName)
		return tmpFileName

	def openCurve(self,fileName):

		if self.instrument == "LATAP":
			desc = ('Time','Rate','Error','Time_Err','Rate_Err2','Exposure')
			ds = 3

		elif self.instrument == "BAT":

			desc = ('Time','Rate','Error','Year','Day','Stat_Err','Sys_Err','Flag',
				'TIMEDEL','TIMEDEL_CODED','TIMEDEL_DITHERED')
			ds = 5

		data = atable.read(fileName, data_start=ds, Reader=atable.NoHeader, guess=False, names=desc)
		return data

	def weighted_avg_and_std(self,values, weights):
		"""                                                                                                                 
		Returns the weighted average and standard deviation.                                                                
		                                                                                                                    
		values, weights -- Numpy ndarrays with the same shape.                                                              
		"""
		average = np.average(values, weights=weights)
		variance = np.dot(weights, (values-average)**2)/weights.sum()  # Fast and numerically precise                       
		return (average, np.sqrt(variance))

	def calcSourceStatistics(self,data):
	    wmean, wstd = self.weighted_avg_and_std(data['Rate'], weights=1./data['Error'])
	    quad_err = np.sqrt(data['Error']**2 + wstd**2)
	    deviations = (data['Rate'] - wmean)/quad_err
	    return {'deviations':deviations,'time':data['Time']}

	def loadSourceList(self,sourceListFile='slist.dat'):

		if(self.instrument == "LATAP"):
		    self.sourceList = atable.read(sourceListFile, names=("name", "ra", "dec"))
		elif(self.instrument == "BAT"):
			tempFile = urlopen('http://swift.gsfc.nasa.gov/results/transients/')
			page = tempFile.read()
			soup = BeautifulSoup(page)
			table = soup.find("table")	
			data = [row.findAll("td") for row in table.findAll("tr")[1:]]
			desc = ('id','name','ra','dec','alt_name','type')
			types =('i4', 'a128', 'f8','f8','a64','a64')
			self.sourceList = Table(names=desc,dtypes=types)
			for row in data:
				cells = [row[0].text,
					str(row[1].find("a")['href']),
					row[2].text,
					row[3].text,
					row[4].text,
					row[5].text]
				self.sourceList.add_row(cells)
		else:
			print "Unknown instrument."
			self.sourceList = -1

		if(self.verbose):
			print "Looking for {} sources.".format(len(self.sourceList))
	



	def loadSourceData(self,dataDir,useWeb=False):

	    if(self.verbose):
		    print "Loading data from {}.".format("the web" if useWeb else dataDir)
	    all_data = [{'name':name, 'data':self.calcSourceStatistics(self.openCurve(self.getCurve(name, directory=dataDir,web=useWeb)))} 
	    	for name in self.sourceList['name']]
	    if(self.verbose):
		    print "Loaded data for {} sources.".format(len(all_data))
	    return all_data

	def findFlares(self,sourceData,days_back,std_above):

	    data_min =  min(np.hstack([source['data']['deviations'] for source in sourceData]))
	    data_max =  max(np.hstack([source['data']['deviations'] for source in sourceData]))
	    data_mean = np.mean(np.hstack([source['data']['deviations'] for source in sourceData]))
	    data_std =  np.std(np.hstack([source['data']['deviations'] for source in sourceData]))


	    if(self.verbose):
		    print "The largest negative deviation is {} cm^-2 s^-1.".format(data_min)
		    print "The largest positive deviation is {} cm^-2 s^-1.".format(data_max)
		    print "The mean deviation is {} cm^-2 s^-1.".format(data_mean)
		    print "The standard deviation of the deviations is {} cm^-2 s^-1.".format(data_std)
		    print "Looking for flares greater than {} cm^-2 s^-1 starting on MJD {}.".format(
		    	data_std*std_above,
		    	self.mjd - days_back)
	    
	    flares = []
	    for source in sourceData:
	        mask = (source['data']['deviations'] > std_above*data_std) & (source['data']['time'] > self.mjd - days_back)
	        if(mask.any()):
	            flares.append( {'name':source['name'], 'time':source['data']['time'][mask], 'deviations':source['data']['deviations'][mask]} )
	    return flares

	def getAssoc(self,sourceNames):

		assocs = []

		if self.instrument == "LATAP":
		
			tmpFileName,response =  urlretrieve("http://fermi.gsfc.nasa.gov/ssc/data/access/lat/2yr_catalog/gll_psc_v08.fit")
			hduFile = pyfits.open(tmpFileName)	

			if(self.verbose):
				print "Loaded 2FGL from FSSC."
				print "Stored the 2FGL fits file as {}.".format(tmpFileName)
			
			for sourceName in sourceNames:
				mask = (hduFile[1].data.field('Source_Name') == re.sub('FGLJ', 'FGL J', sourceName))
				sourceInfo = (hduFile[1].data)[mask]
				assocs.append({'name':sourceName,'assoc':sourceInfo.field('ASSOC1')[0],'class':sourceInfo.field('CLASS1')[0]})
		
		elif self.instrument == "BAT":

			for sourceName in sourceNames:
				mask = self.sourceList['name'] == sourceName
				sourceInfo = self.sourceList[mask]
				name = re.sub('weak/','',self.sourceList[mask]['name'][0])
				assocs.append({'name':name,'assoc':name,'class':self.sourceList[mask]['type'][0]})

		else:
			print "Unknown instrument."
			return -1

		return assocs
	
	def getLC(self,sourceName,short=False):

		if self.instrument == "LATAP":

			if short:
				URL = "ap_lcs/lightcurve_" + re.sub('\+', 'p', sourceName) + ".png"
			else:
				URL = self.baseURL['LATAP'] + "lightcurve_" + re.sub('\+', 'p', sourceName) + ".png"

		elif self.instrument == "BAT":

			URL = self.baseURL['BAT'] + sourceName + ".png"

		else:

			print "Unknown Instrument"
			URL = -1

		return URL

	def simbadByName(self,sourceName):

		URL = "http://simbad.u-strasbg.fr/simbad/sim-id?Ident="

		if self.instrument == "LATAP":
			sourceName = re.sub('\+', '%2B', sourceName)	
			sourceName = re.sub(' ', '+', sourceName)
		return URL + sourceName

	def simbadByCoord(self,ra,dec):
		
		URLprefix = "http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord="
		URLsuffix = "&submit=submit%20query&Radius.unit=arcmin&CooEqui=2000&CooFrame=FK5&Radius=10"
		return "{}{}%20{}{}".format(URLprefix,ra,dec,URLsuffix)

	def dataByName(self,sourceName):
		
		if self.instrument == "LATAP":
			URL =  self.baseURL['LATAP'] + "lc_" + re.sub('\+', 'p', sourceName) + ".dmp1.out"

		elif self.instrument == "BAT":

			URL = self.baseURL['BAT'] + sourceName

		else:

			print "Unknown Instrument"
			URL = -1

		return URL

	def printFlares(self,flares,html,htmlfile,days_back):

		sL = self.sourceList
		assocs = self.getAssoc([flare['name'] for flare in flares])
		
		html_string = '<p class="center"><b>Looking for flares from MJD {} to MJD {}.</b></p>'.format(self.mjd-days_back,self.mjd)
		html_string += '<table class="styled-table" width="100%">\n'
		count = 0
		for flare in flares:
			for assoc in assocs:
				matched = assoc['name'] == re.sub('weak/','',flare['name'])
				this_assoc = assoc['assoc'] if matched else 0
				if matched:
					break
		
			print 
			print flare['name']
			print "\tAssociation: {}".format(this_assoc)
			print "\tLocation: {},{}".format(sL[sL['name'] == flare['name']]['ra'][0],sL[sL['name'] == flare['name']]['dec'][0])
			print "\tFlare Level(s): {}".format(flare['deviations'])
			print "\tFlare Date(s): {}".format(flare['time'])
			print "\tLight Curve: {}".format(self.getLC(flare['name'],False))

			if(np.mod(count,3) == 0):
				html_string += '<tr>\n'
			html_string += '<td class="center">\n'
			html_string += '<b>{}</b><br>\n'.format(flare['name'])
			html_string += '<a href="{}">({},{})</a><br>\n'.format(
				self.simbadByCoord(sL[sL['name'] == flare['name']]['ra'][0],sL[sL['name'] == flare['name']]['dec'][0]),
				sL[sL['name'] == flare['name']]['ra'][0],
				sL[sL['name'] == flare['name']]['dec'][0])
			fileName = self.getLC(flare['name'],False)
			html_string += '<a href="{}"><img src="{}" width="200" height="155" alt="{}"></a><br>\n'.format(fileName,fileName,fileName)
			html_string += 'Data: <a href="{}">{}</a><br/>\n'.format(self.dataByName(flare['name']),flare['name'])
			html_string += 'Assoc: <a href="{}">{}</a>'.format(self.simbadByName(this_assoc),this_assoc)

			#html_string += <a href="ap_lcs/lc_2FGLJ0000.9-0748.dmp1.out">Data</a>
			html_string += '</td>\n'
			if(np.mod(count+1,3) == 0):
				html_string += '</tr>\n'
			
			count = count + 1
		print ""
		#Finish off the table if you're not on a cell divisible by 3.
		if(np.mod(count,3) != 0):
			for n in range(3 - np.mod(count,3)):
				html_string += '<td class="center">Intentionally Left Blank.</td>\n'
			html_string += "</tr>\n"
		html_string += '</table>\n'
		html_string += '<p class="center">Generated using version {}.</p>'.format(__version__)
		if(html):
			if(self.verbose):
				print "Writing html file {}.".format(htmlfile)
			f = open(htmlfile, 'w')
			f.write(html_string)
			f.close()

	def FlareFinder(self, days_back, std_above, sourceListFile, dataDir, fromWeb, html, htmlfile):

	    if(self.verbose):
		    print "This is version {}.".format(__version__)

	    self.loadSourceList(sourceListFile)
	    data = self.loadSourceData(dataDir,fromWeb)
	    flares = self.findFlares(data,days_back,std_above)
	    self.printFlares(flares,html,htmlfile,days_back)

def cli():

	helpString = "Looks at the adaptive binned lightcurves for flares.\
                  It first imports all of the lightcurve data, figures\
                  out the mean and standard deviation of the data and \
                  then prints out flares above a certain standard deviation\
                  within a certain time period."

	import argparse

	parser=argparse.ArgumentParser(description=helpString)

	parser.add_argument("--instrument", type=str, default="LATAP", help="LATAP or BAT")
	parser.add_argument("--mjd", type=int, default=0, help="The starting MJD (if 0, will use current).")
	parser.add_argument("--days_back", type=int, default=28, help="The lookback time in days (default = 28).")
	parser.add_argument("--threshold", type=float, default=4, help="The flare threshold in standard deviations (default = 4).")
	parser.add_argument("--sourceList", default="slist.dat", help="Text file containing the names of the sources.")
	parser.add_argument("--dataDir", default="lcs", help="Directory containing the lightcurve files (default = lcs).")
	parser.add_argument("--fromWeb", default=False, help="Read the lightcurve data from the FSSC website or not.")
	parser.add_argument("--html", default=False, help="Print an html file.")
	parser.add_argument("--htmlfile", default='flares.html', help="Output html file (default = flares.html).")
	parser.add_argument("--verbose", default=False, help="Turn on or off the debug statements.")


	args = parser.parse_args()

	ff = FlareFinder(verbose=args.verbose, 
		instrument=args.instrument,
		mjd = args.mjd)

	ff.FlareFinder(args.days_back, 
		args.threshold, args.sourceList, 
		args.dataDir, args.fromWeb, 
		args.html, args.htmlfile)

if __name__ == '__main__': cli()
