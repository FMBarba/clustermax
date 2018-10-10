

*****************************************************************************
* Clustermax example: 400 random points between Illinois and Indiana, USA.
*****************************************************************************


* Generate 400 random points on a 150 x 450 km area
	set obs 400
	set seed 1
	gen lat = 40  + runiform()
	gen lon = -90 + runiform()*4

* Set minimum cluster size to one; minimum between-cluster distance to 25km; maximum within-cluster distance to 15km:	
	clustermax lat lon, seed(1) gen(cluster) w(15) b(25) n(1)

