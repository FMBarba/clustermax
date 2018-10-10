{smcl}
{title:Title}

{phang}
{bf:clustermax} {hline 2} Create a maximum number of clusters from given point coordinates


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:clustermax}
{it:lat lon}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(varname)}}the new cluster id variable{p_end}
{synopt:{opt n(#)}}minimum number of points per cluster{p_end}
{synopt:{opt w:ithin(#)}}maximum distance between any two points within the same cluster{p_end}
{synopt:{opt b:etween(#)}}minimum distance between any two points from different clusters{p_end}
{synopt:{opt seed(#)}}seed value{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:clustermax} groups point coordinates into a maximum feasible number of clusters given user-specified cluster characteristics. 

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt n} the minimum number of points within cluster.

{phang}
{opt w:ithin} points that form a cluster fit into a circle of diameter {it: within}.{it: within} <= {it: between} is required to uniquely assign points to clusters. {it:within}=0 is allowed for {it:n}=1.

{phang}
{opt b:etween} the minimum distance between any two points from different clusters. 

{phang}
{opt seed} ensures the same points are selected each time when the choice between two points irrelevant to the maximization problem. 

{marker example}{...}
{title:Examples}

{pstd}Generate 100 random points on a 110 x 340 km area in Illinois, USA.{p_end}
	{cmd:. set obs 100}
	{cmd:. set seed 1}
	{cmd:. gen lat = 40  + runiform()}
	{cmd:. gen lon = -90 + runiform()*4}

{pstd}Set cluster size to one; minimum between-cluster distance 20km:{p_end}
	{cmd:. clustermax lat lon, seed(1) gen(cluster) b(20) n(1)}

{pstd}Set minimum cluster size to one; minimum between-cluster distance 20km; maximum within-cluster distance 20km:{p_end}
	{cmd:. clustermax lat lon, seed(1) gen(cluster) w(20) b(20) n(1)}

{pstd}Set minimum cluster size to three; minimum between-cluster distance 25km; maximum within-cluster distance 15km:{p_end}
	{cmd:. clustermax lat lon, seed(1) gen(cluster) w(15) b(25) n(3)}


{title:Required packages}

{pstd}{cmd:clustermax} requires the packages {it: geodist, tuples, distinct} and {it: ftools}, available on SSC.{p_end}

{marker Aknowledgements}{...}
{title:Aknowledgements}

{pstd}Jef L Leroy provided useful comments on this command.{p_end}

{marker Author}{...}
{title:Author}

{pstd}Francisco M Barba, International Food Policy Research Institute, Washington DC, USA{p_end}

