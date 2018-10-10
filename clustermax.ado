*********************
* clustermax.ado
* version 9
* Francisco Barba
* October 8, 2018
*********************

cap program drop clustermax
program clustermax, rclass

	syntax varlist , seed(integer) GENerate(name) n(integer) [Within(real 0)]  Between(real) 
	qui {
	qui {
		*************************************************************************************************************************
		* PART 1: Find all combinations of points with MUTUAL maximum distance of `within'
		*************************************************************************************************************************
		
		cap gen `generate' 
		if _rc==110 {
			di as err "variable `generate' already defined"
			exit 110
		}	
		if _rc==0 {
			cap drop `generate' 
		}		
		
		
		set seed `seed'
		
		tempvar sorttemp
		gen double `sorttemp' = runiform()	
		sort `sorttemp'
		tempvar id
		qui  gen `id'=_n	
		tempvar id_me
		qui  gen `id_me'=`id'		
		qui tempfile original0
		qui save `original0' 	
		
		local clustervarname `generate' 
		token `varlist'
		local latitude  `1'
		local longitude `2'
		* Define errors
		if `n'>1 & "`within'"=="" {
			di as errr "Error: must specift within if n>1"
			continue, break
		}
		if `within'>`between'{
			di as err "Error: between() must be greater or equal witin_dist()."
			continue, break
		}
		if (`between'<=0) { //| 
			di as err "Error: between()  must be positive."
			continue, break
		}	
		if `within'<0  { 
			di as err "Error: within()  must be greater/equal zero."
			continue, break
		}	
		if `within'==0 & `n'>1  { 
			di as err "Error: within=0 only allowed for n=1."
			continue, break
		}			
		keep `id' `latitude' `longitude'

		tempfile original
		save `original', replace
		
		tempvar dup3
		qui duplicates tag `latitude' `longitude', gen(`dup3')
		qui  fsum `dup3'
		if `r(max)'!=0 {
			di as err "Error: duplicate lat/lon coordinates found. lat/lon must be unique."
			exit
		}		
		cap drop `dup3'
		tempvar dup
		qui duplicates tag `latitude' , gen(`dup')
		}	// ends qui
		qui  fsum `dup'
		if `r(max)'!=0 {
			di as text "Warning: duplicate lat coordinate found."
		}
		cap drop `dup'
		tempvar dup2
		qui duplicates tag `longitude' , gen(`dup2')
		qui fsum `dup2'
		if `r(max)'!=0 {
			di as text "Warning: duplicate lon coordinate found."
		}
		cap qui drop `dup2'
		qui count if `latitude'==.
		if `r(N)'!=0 {
			di as  "Warning: missing lat coordinate found."
		}
		qui count if `longitude'==.
		if `r(N)'!=0 {
			di as  "Warning: missing lon coordinate found."
		}
		qui {
			tempvar latitude0 
			tempvar longitude0
			tempvar id0
			* Skip the loop that identifies neighbors if `within'=0 (allowed when n=1) 
			if `within'!=0 {	
				rename `latitude'  `latitude0' 
				rename `longitude' `longitude0'
				rename `id' `id0'

				cross using `original' 

				tempvar d
				geodist `latitude0' `longitude0' `latitude' `longitude', gen(`d')
				
				* Get neighbors (for later)
				
					tempvar x2
					bys `id0': gen `x2'=`d'<=`within' if `d'!=0
					tempvar points_close2
					bys `id0': egen `points_close2'=sum(`x2')					
					* save number of points in potential conflict area
					tempvar y
					bys `id0': gen `y' =`d'<=`between' if `d'!=0 // &d>`within'
					tempvar points_within_ring
					bys `id0': egen `points_within_ring'=sum(`y')				
					drop `x2'
					
				preserve	
					keep `id0' `points_close2' `points_within_ring'
					duplicates drop `id0', force
					tempfile temp_points_close
					save `temp_points_close'
				restore	
					
				drop `points_close'  `points_within_ring'
				
				* at least `n'-1 points in inner cluster 
				tempvar x
				bys `id0': gen `x'=`d'<=`within' if `d'!=0
				tempvar points_close
				bys `id0': egen `points_close'=sum(`x')			
				drop `d' `x'	
			
				* This while-loop keeps dropping observations until every observation has a minimum number of `n' neighbor points. 
				//It needs to be done iteratively:
				// example: p1 is (only) close to p2 and p3, but p2 is not close to p3
				// hence, after dropping p2 or p3 (which only have 1 as a neighbor) p1 is only left with one neighbor.
				// hence, the dropping of iterations has to be iterated again, until EACH point has at least `n' neighbors			
				fsum `points_close'
				if `r(max)'!=0 {
					local min_neigbor_points=-1
					while `min_neigbor_points'<`n'-1   {
						keep if `points_close'>=`n'-1 

						count 
						if `r(N)' ==  0  {
							di as err "Error: There don't exist `n' points with a max. pairwise distance of `within' in the current data."
							exit 111
						}	
						keep `id0' `latitude0' `longitude0'
						qui duplicates drop `id0', force

						tempfile temp1
						qui save `temp1'
						
						rename  `latitude0' `latitude'
						rename  `longitude0' `longitude'
						rename  `id0' `id'
						cross using `temp1'
						
						cap drop d
						tempvar d
						geodist `latitude0' `longitude0' `latitude' `longitude', gen(`d')
						
						tempvar x
						bys `id0': gen `x'=`d'<=`within' if `d'!=0
						tempvar points_close
						bys `id0': egen `points_close'=sum(`x')
						fsum `points_close'
						local min_neigbor_points=`r(min)'
					}
					keep if `x'==1 | `x'==.
				}
			}	
			* Match potential points for clusters to each other (only if at least 2 points within cluster)
			* seed point has at least `n'-1 neighbors 		
			if `n' > 1 {	
				****************
				* IMPORTANT SORT! Establishes reshape order
				sort `id0' `id' 
				****************
				keep `id0' `id'
				tempvar nnn
				bys `id0': gen `nnn'=_n
				* This reshape gives the neighbor IDs in a row next to the source ID
				reshape wide `id', i(`id0') j(`nnn')
				
				local cols=`r(k)'-1 
				local c string(`id'1)
				forval i= 2 / `cols' {
					* This appending macro loop generates a string of original ID plus all neighbor IDs, separated by two blanks
					local c `c' + " " + " " + string(`id'`i')
				}
				local c `" "  " + " " +  `c' + " " + " " "'
				tempvar comb
				gen `comb'=  `c'
				replace `comb'=subinstr(`comb', ".", "", .)

				cap drop `nnn'
				* Row counter
				tempvar nnn
				gen `nnn'=_n
				* This loop starts the matching process of dyads, triads... if n is 2, 3... correspondingly
				* Start with first row tuple
				count
				forval src_row=1/`r(N)' { // foreach seed (original point)
					* get the source seed value in `src_row'
					local src_seed=`id0' in `src_row'	
					* get neighbor combination string (including seed)  of `src_row'	
					local src_comb_`src_seed'=`comb' in `src_row'
					*find the column number of the seed within the string: (cond() value of the tuple) 
					local tuplelength: word count `src_comb_`src_seed''
					forval j=1/`tuplelength' { 
						local column 
						local num`j': word `j' of `src_comb_`src_seed''
						if `num`j''==`src_seed' {
							local column`src_seed'=`j'
						}
					}
					* run all tuples
					cap tuples `src_comb_`src_seed'', display min(`n') max(`n') cond(`column`src_seed'')
					if _rc==3300 {
						di as err "Error: Mata out of memory. Too many neighbors to evaluate. Consider decreasing the within-cluster distance or n."
						exit 3300
					}
					* Get the number of tuples of `src_row'
					local count=" "
					local num=1
					while "`count'"!="" {
						if "`tuple`num''"=="" {
							local count=""
						}
						local src_tuple`num' `tuple`num''
						local tuple`num'
						local ++num
					}
					local num_tuples_`src_seed'=`num'-2
					
					* Get neighbor elements of the each tuple of `src_seed'
					forval src_tuple= 1/`num_tuples_`src_seed'' { // for all tuples of `src_seed'
						local src_`src_seed'_`src_tuple'="`src_tuple`src_tuple''"
						local nbrs_`src_seed'_`src_tuple'= subinstr(" `src_tuple`src_tuple'' ", " `src_seed' ", " ", .) // the n-1 neighbor seeds 
					
						local nbrseed=1
						* Get the seed `nbrseed' of `src_tuple' in `src_seed'
						local n_within_cluster_min1=`n'-1
						forval nbrseed=1/`n_within_cluster_min1'{
							local nbr_seed : word `nbrseed' of `nbrs_`src_seed'_`src_tuple''
							* Get row of neighbor seed 
							
							cap drop `x'
							tempvar x
							qui gen `x'=`nnn' if `id0'==`nbr_seed'
							qui fsum `x'
							local nbr_row=`r(max)'
							drop `x'
							
							* get neighbor combination string (including seed)  of `nbr_row'	
							local nbr_comb_`src_seed'_`src_tuple'_`nbr_seed' =`comb' in `nbr_row'
							*find the column number of the seed within the string: (cond() value of the tuple) 
							local tuplelength: word count `nbr_comb_`src_seed'_`src_tuple'_`nbr_seed''
							forval j=1/`tuplelength' { 
								local column 
								local num`j': word `j' of `nbr_comb_`src_seed'_`src_tuple'_`nbr_seed''
								if `num`j''==`nbr_seed' {
									local column`nbr_seed'=`j'
								}
							}
							* Tuples of neighbor seed
							cap tuples  `nbr_comb_`src_seed'_`src_tuple'_`nbr_seed'', display min(`n') max(`n') cond(`column`nbr_seed'') 
							if _rc==3300 {
								di as err "Error: Too many neighbors to evaluate. Consider decreasing the within-cluster distance"
								exit 3300
							}
							* Get the number of tuples of `nbr_row'
							local count2=" "
							local num2=1
							while "`count2'"!="" {
								if "`tuple`num2''"=="" {
									local count2=""
								}
								local nbr_tuple`num2' `tuple`num2''
								local tuple`num2'
								local ++num2
							}
							local num_tuples2=`num2'-1	
		}
							* If the tuple of that neighbor seed is the same than the one from the originating source seed
							* Generate a match and save the IDs into a local of matched tuples
							forval i = 1/`num_tuples2' {
								if "`src_`src_seed'_`src_tuple''"=="`nbr_tuple`i''" {
									if `src_row'==1 & `src_tuple'==1 {
										noisily di as text "matching neighbors..."
									}
									local matched_tuples `"`matched_tuples' "`src_`src_seed'_`src_tuple''""' 
									local tuple`num2'
								}	
							}		
						}
					}
				
				* After the matching loop ran for each row, save the string of matched tuples into a dataset
				local n_used_tuples: word count "`matched_tuples'"
				preserve
					clear all
					set obs  `n_used_tuples'
					cap drop `clus'
					tempvar clus
					gen `clus'=""
					* Select one tuple 
					forval xy=1/`n_used_tuples'{
						local t: word `xy' of `matched_tuples'
						replace `clus' ="`t'" in `xy'
					}
					* Save the tuple number as cluster id into dataset
					fegen `clustervarname'=group(`clus')	
					cap drop `NNN'
					tempvar NNN
					bys `clustervarname': gen `NNN'=_N
					drop if `clustervarname'==. // check why missings created
					* keep only tuples that matched at least `n' times, e.g. if `n'=3 :  
					// 1 neighor of 2 and 3; 2 neighor of 1 and 3; 3 neighor of 1 and 2; : 3 matches				
					keep if `NNN'>=`n'

					collapse (firstnm) `clus' (mean) `NNN' , by(`clustervarname')
					 
					  * prepare reshape
					 tempvar clus_iid
					 split `clus', parse(" ") generate(`clus_iid')
					 destring `clus_iid'*, replace
					 
					 drop `clus' `NNN'
					 tempvar new
					 * reshape to long (one id per row, several rows per cluster ID)				 
					 reshape long `clus_iid', i(`clustervarname') j(`new')
					 rename `clus_iid' `id0'
					 
					 tempfile temp2
					 save `temp2', replace
				restore	
				merge 1:m `id0' using `temp2', nogen

				keep `id0' `clustervarname'
				rename `id0' `id'

				merge m:1 `id' using  `original'
				keep if _merge==3
				drop _merge

			}
		}
		* Special case: if only one cluster feasible in total, skip cluster conflict resolution section
		if `n' != 1 {
			distinct `clustervarname' 
		}
		if `n' == 1 {
			distinct `id' 
		}		
		local distinct_clusters `r(ndistinct)'
		if `distinct_clusters'!= 1 {
		
			*************************************************************************************************************************
			* PART 2: Resolve between-cluster conflicts
			*************************************************************************************************************************
			if `n' == 1 & `within'!=0 {
				duplicates drop `id', force
				drop `latitude0' `longitude0' 
				drop `id0'
			}
			tempfile temp3
			qui save `temp3'	
			if `n' != 1 {
				rename `clustervarname' `clustervarname'0	
			}	
			rename  `latitude' `latitude0'
			rename  `longitude' `longitude0'
			rename  `id' `id0'
			
			cross using `temp3'
			cap drop d
			tempvar d
			geodist `latitude0' `longitude0' `latitude' `longitude', gen(`d')		
			
			drop if `d'==0
			if `n' != 1 {
				drop if `clustervarname'0 == `clustervarname' 
			}
			if `n' == 1 {
				drop if `id0' == `id' 
			}			
			if `n' != 1 {
				drop if `clustervarname' ==. | `clustervarname'0 ==.
			}			
			if `n' == 1 {
				drop if `id' ==. | `id0' ==.
			}
			count if `d'<`between' 
			local tot_conflicts `r(N)'
			
			* indicator if any conflict between points of different clusters
			tempvar conflict_distance
			gen `conflict_distance'=`d'<`between'

			* Save dataset with id sorting from seed 
			gen `id_me' = `id0'
			merge m:1 `id_me' using `original0', keepusing( `sorttemp')
			cap drop `id_me' 
			if `n' != 1 {
				sort  `clustervarname'0 , stable
			}
			tempvar sortorder
			gen double `sortorder' = runiform()	
			tempvar sort_id
			bys `id0': egen double `sort_id'=min(`sorttemp') 	

			preserve
				qui duplicates drop `id0', force
				keep `id0' `sort_id'
				//rename `id' `id0'
				tempfile temp_sort_id
				qui save `temp_sort_id'
			restore
			
			* Save dataset with cluster sorting from seed 
			if `n'!=1 {
				tempvar min_sortorder
				bys `clustervarname'0: egen double `min_sortorder'=min(`sortorder') 	
				tempvar sort_clustervarname
				fegen `sort_clustervarname'=group(`min_sortorder')
				
				preserve
					qui duplicates drop `clustervarname'0, force
					keep `clustervarname'0 `sort_clustervarname'
					//rename `clustervarname' `clustervarname'0
					cap drop _merge
					tempfile temp_sort_clustervarname
					qui save `temp_sort_clustervarname'
				restore
			}					
			**********************************************
			* loop for conflict resolution: initialization 	
			if `n' != 1 {
				cap drop `a'
				tempvar a
				bys    `clustervarname'0   `clustervarname' `conflict_distance': gen `a'=_n if `conflict_distance'==1 & `clustervarname'!=.
				merge m:1 `id0' using `temp_points_close', nogen  // sort by number of points close
				tempvar cluster_points_close
				bys `clustervarname'0: egen `cluster_points_close'=sum(`points_close2')
				tempvar rank_cluster_points_close			
				fegen `rank_cluster_points_close' =group(`cluster_points_close')
				fsum `rank_cluster_points_close'
				tempvar rev_rank_cluster_points_close
				gen `rev_rank_cluster_points_close'=`r(max)'+`r(min)'-`rank_cluster_points_close'
				tempvar rank_points_within_ring
				fegen `rank_points_within_ring' =group(`points_within_ring')
				fsum `rank_points_within_ring'
				tempvar rev_rank_points_within_ring
				gen `rev_rank_points_within_ring'=`r(max)'+`r(min)'-`rank_points_within_ring'
				merge m:1 `clustervarname'0 using `temp_sort_clustervarname' , nogen
			}								
			if `n' == 1 {
				cap drop `a'
				tempvar a
				bys    `id0'   `id' `conflict_distance': gen `a'=_n if `conflict_distance'==1 & `id'!=.		
				if `within'!=0 {
					merge m:1 `id0' using `temp_points_close', nogen  // sort by number of points close
					tempvar cluster_points_close
					bys `id0': egen `cluster_points_close'=sum(`points_close2')
					tempvar rank_cluster_points_close			
					fegen `rank_cluster_points_close' =group(`cluster_points_close')
					fsum `rank_cluster_points_close'
					tempvar rev_rank_cluster_points_close
					gen `rev_rank_cluster_points_close'=`r(max)'+`r(min)'-`rank_cluster_points_close'
					tempvar rank_points_within_ring
					fegen `rank_points_within_ring' =group(`points_within_ring')
					fsum `rank_points_within_ring'
					tempvar rev_rank_points_within_ring
					gen `rev_rank_points_within_ring'=`r(max)'+`r(min)'-`rank_points_within_ring'
				}
				merge m:1 `id0' using `temp_sort_id', nogen		
			}	
			tempvar c
			gen `c'=`a'==1	
			if `n' != 1 {
				tempvar cluster_conflicts
				bys     `clustervarname'0  : egen `cluster_conflicts'= sum(`c')
			}	
			if `n' == 1 {
				tempvar cluster_conflicts
				bys    `id0'  : egen `cluster_conflicts'= sum(`c')
			}	
			fsum `cluster_conflicts' if `cluster_conflicts'>0
			if `r(N)'== 0  {
				local no_conflicts 1
			}
			if `r(N)'!= 0  {  // if there is any conflict: (otherwise go directly to end)
				if `n' != 1 {
					local max_temp=`r(max)'
					local min_temp=`r(min)'
					tempvar neighbor_to_min

					gen `neighbor_to_min'=0
					fsum `cluster_conflicts'   if `cluster_conflicts'>0
					local min_greater_one=`r(min)'
					tempvar n_neighbors_min // number  of neighors that have `cluster_conflicts'==`r(min)'
					tempvar xx 
					bys `clustervarname'0 `clustervarname' `conflict_distance': gen `xx' =_n if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1
					tempvar xx2
					gen `xx2' =`xx' ==1
					bys     `clustervarname'0  : egen `n_neighbors_min'= sum(`xx2') 
					flevelsof `clustervarname' if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1, local(levels)
					foreach v in `levels' {
						replace `neighbor_to_min'=1 if `clustervarname'0==`v'
					}	
					tempvar sortorder_cluster
					fegen `sortorder_cluster' = group(  `neighbor_to_min' `n_neighbors_min' `sort_clustervarname') //if `e'==1	
				}
				if `n' == 1 {
					local max_temp=`r(max)'
					local min_temp=`r(min)'
					tempvar neighbor_to_min
					gen `neighbor_to_min'=0
					fsum `cluster_conflicts'   if `cluster_conflicts'>0
					local min_greater_one=`r(min)'
					tempvar n_neighbors_min // number  of neighors that have `cluster_conflicts'==`r(min)'
					tempvar xx 
					bys `id0' `id' `conflict_distance': gen `xx' =_n if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1
					tempvar xx2
					gen `xx2' =`xx' ==1
					bys     `id0'  : egen `n_neighbors_min'= sum(`xx2') 
					flevelsof `id' if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1, local(levels)
					foreach v in `levels' {
						replace `neighbor_to_min'=1 if `id0'==`v'
					}	
					tempvar sortorder_cluster
					fegen `sortorder_cluster' = group( `neighbor_to_min'  `n_neighbors_min' `sort_id') //if `e'==1
				}	
				
				if `n' != 1 {
					tempvar max_conflict_distance
					bys `clustervarname'0 `clustervarname': egen `max_conflict_distance'=max(`conflict_distance')
					duplicates drop `clustervarname'0 `clustervarname', force
					replace `conflict_distance'=`max_conflict_distance'
					drop `max_conflict_distance' `id' `d' `latitude' `longitude'

					fsum `sortorder_cluster'
					if `r(N)'!=0 {
						flevelsof `clustervarname'0 if `sortorder_cluster'==`r(max)' , local(levels)
						recode `conflict_distance'(1=0) if `clustervarname'==`levels'
						drop if `clustervarname'0 == `levels' 
					}	
				}
				if `n' == 1 {
					tempvar max_conflict_distance
					bys `id0' `id': egen `max_conflict_distance'=max(`conflict_distance')
					duplicates drop `id0' `id', force
					replace `conflict_distance'=`max_conflict_distance'
					drop `max_conflict_distance' `d' `latitude' `longitude' 

					fsum `sortorder_cluster'
					flevelsof `id0' if `sortorder_cluster'==`r(max)' , local(levels)
					
					recode `conflict_distance'(1=0) if `id'==`levels'
					drop if `id0' == `levels' 
				}
			}

		**********************************************************************************************************
		* loop for conflict resolution: loop through all conflicts
		
			fsum `cluster_conflicts'
			if `r(N)'!=0 {  // if there is any conflict: (otherwise go directly to end)
				local max_conflicts=`r(max)'
				while `max_conflicts'>0 {			
					qui {
						if `n' != 1 {
							cap drop   `cluster_conflicts'
							cap drop `a'
							tempvar a	
							bys    `clustervarname'0   `clustervarname' `conflict_distance': gen `a'=_n if `conflict_distance'==1 & `clustervarname'!=.
							tempvar c
							gen `c'=`a'==1	
							tempvar cluster_conflicts
							bys    `clustervarname'0  : egen `cluster_conflicts'= sum(`c')
						}
						if `n' == 1 {
							cap drop   `cluster_conflicts'
							cap drop `a'
							tempvar a	
							bys    `id0'   `id' `conflict_distance': gen `a'=_n if `conflict_distance'==1 
							tempvar c
							gen `c'=`a'==1	
							tempvar cluster_conflicts
							bys    `id0'  : egen `cluster_conflicts'= sum(`c')
						}				
						fsum `cluster_conflicts' if `cluster_conflicts'>0
						if `r(N)'!= 0  {  //	
							if `n' != 1 {
								drop `neighbor_to_min' `n_neighbors_min' 
								tempvar neighbor_to_min
								gen `neighbor_to_min'=0
								fsum `cluster_conflicts'   if `cluster_conflicts'>0
								local min_greater_one=`r(min)'
								tempvar n_neighbors_min // number  of neighors that have `cluster_conflicts'==`r(min)'
								tempvar xx 
								bys `clustervarname'0 `clustervarname' `conflict_distance': gen `xx' =_n if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1
								tempvar xx2
								gen `xx2' =`xx' ==1
								bys     `clustervarname'0  : egen `n_neighbors_min'= sum(`xx2') 
								flevelsof `clustervarname' if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1  , local(levels)
								foreach v in `levels' {
									replace `neighbor_to_min'=1 if `clustervarname'0==`v'
								}				
								drop `sortorder_cluster' 
								tempvar sortorder_cluster
								fegen `sortorder_cluster' = group(`neighbor_to_min' `n_neighbors_min' `sort_clustervarname') //if `e'==1	
								fsum `sortorder_cluster'
								if `r(N)'!=0 {
									flevelsof `clustervarname'0 if `sortorder_cluster'==`r(max)'  , local(levels)
									local clustervarname_todrop= `levels'
									recode `conflict_distance'(1=0) if `clustervarname'==`clustervarname_todrop'
									drop if `clustervarname'0 == `clustervarname_todrop'
								}	
							}
							if `n' == 1 {
								drop `neighbor_to_min' `n_neighbors_min' 
								tempvar neighbor_to_min
								gen `neighbor_to_min'=0
								fsum `cluster_conflicts'   if `cluster_conflicts'>0
								local min_greater_one=`r(min)'
								tempvar n_neighbors_min // number  of neighors that have `cluster_conflicts'==`r(min)'
								tempvar xx 
								bys `id0' `id' `conflict_distance': gen `xx' =_n if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1
								tempvar xx2
								gen `xx2' =`xx' ==1
								bys     `id0'  : egen `n_neighbors_min'= sum(`xx2') 
								flevelsof `id' if `cluster_conflicts'==`min_greater_one' & `conflict_distance'==1  , local(levels)
								foreach v in `levels' {
									replace `neighbor_to_min'=1 if `id0'==`v'
								}				
								drop `sortorder_cluster' 
								fegen `sortorder_cluster' = group( `neighbor_to_min'  `n_neighbors_min' `sort_id') //if `e'==1	
								fsum `sortorder_cluster'
								flevelsof `id0' if `sortorder_cluster'==`r(max)'  , local(levels)
								local clustervarname_todrop= `levels'
								recode `conflict_distance'(1=0) if `id'==`clustervarname_todrop'
								drop if `id0' == `clustervarname_todrop'
							}					
					} // ends qui
							fsum `cluster_conflicts'
							local max_conflicts_print `max_conflicts' 
							if `r(N)'!= 0  {  
								if `r(max)'==`max_conflicts' {
									//local dots `dots'.
									//local max_conflicts_print "resolving conflicts:`max_conflicts_print'  `dots'"
									//di as input "`max_conflicts_print'" _continue
									//di as input "`dots'" _continue
									//di as err  "max_conflicts=`max_conflicts_print' `dots'" _continue
								}
								*
								if `r(max)'!=`max_conflicts' {
									local dots
									noisily di as text "resolving conflicts:`r(max)'"
								}		
								local max_conflicts=`r(max)'
							}
						}
						
					if `r(N)'== 0  {  
						local max_conflicts 0
						di as input "resolving conflicts:0"
					}
				}	
			}	
			if _rc==2000 {
				di as err "Error: No observations." 
				di as err "Likely cause is a between-cluster distance which is too large."
				exit 2000
			}
			if `n' != 1 {	
				cap duplicates drop `clustervarname'0 `id0' , force 
				keep `clustervarname'0 `id0' 
			}
			if `n' == 1 {	
				cap duplicates drop `id0' , force 
				keep `id0' 
			}
			rename `id0' `id' 
			merge 1:1 `id' using `original0'
			keep if _merge==3
			drop _merge
			if `n' != 1 {	
				keep `clustervarname'0 `id' `latitude' `longitude' 
			}
			if `n' == 1 {	
				keep `id' `latitude' `longitude' 
			}	
		*************************************************************************************************************************
		* Part 3: After maximum number of clusters was chosen, keep adding points that respect constraints to clusters 	
		*************************************************************************************************************************	
			if `within'!=0 {
				qui {
					//drop if `latitude0'==. | `longitude0'==. 
					if `n' == 1 {
						gen `clustervarname'=_n
					}
					if `n' != 1 {	
						rename  `clustervarname'0 `clustervarname'
					}
					merge 1:1  `latitude' `longitude' using `original0', nogen
					
					* Find all remaining points within `within' distance
					tempfile temp4
					save `temp4'
					
					rename  `id'	`id0'
					rename  `clustervarname' `clustervarname'0
					rename	`longitude' `longitude0'
					rename  `latitude' `latitude0'	
					
					cross using `temp4'
					
					tempvar d
					geodist `latitude0'	 `longitude0' `latitude' `longitude'  , gen(`d')
					drop if `id0'==`id'
					
					* Distance from remaining points to old points in clusters (new points in: cluster_id, old points in:cluster_id0)
					tempvar within_distance
					bys `clustervarname'0: gen `within_distance'=1 if `d'<=`within' & `d'!=0 & `clustervarname'0!=. & `clustervarname'==.
					
					tempvar max_within_dist
					bys `clustervarname'0 `id': egen `max_within_dist'=max(`within_distance')
					tempvar max_d_within
					bys `clustervarname'0 `id': egen `max_d_within'=max(`d') if `max_within_dist'==1
					tempvar new_point
					gen `new_point'=`max_d_within'<=`within' // flag for potential new point
					replace `clustervarname' = `clustervarname'0 if `new_point'==1	
				
					keep if `max_d_within'<=`within' | `clustervarname'!=. // keep potential new points and old assigned points

					keep `id' `clustervarname' 	`latitude' `longitude' `new_point'
					qui duplicates drop `id', force
					
					* First drop points that are in conflict with points in existing clusters (would conflict with cluster maximization)
					tempfile temp5
					save `temp5'
						
					rename  `id'	`id0'
					rename  `clustervarname' `clustervarname'0
					rename	`longitude' `longitude0'
					rename  `latitude' `latitude0'	
					rename  `new_point' `new_point'0
					cross using `temp5'
					
					tempvar d
					geodist `latitude0' `longitude0'  `latitude' `longitude' , gen(`d')	
					drop if `id0'==`id'
					
					tempvar conflict
					gen `conflict'=1 if `d'< `between' & `clustervarname'0 !=`clustervarname' & `new_point'0==0 & `new_point'==1
					
					tempvar max_conflict
					bys `id': egen `max_conflict'=max(`conflict')
					tempvar max_new_point
					bys `id': egen `max_new_point'=max(`new_point')
					drop if `max_conflict'==1 &  `max_new_point'==1	
					
					* Second, sequentially drop new points that have most conflicts with other new points
					keep `id' `clustervarname'   `latitude' `longitude' `new_point'
					qui duplicates drop `id', force	
					tempfile temp6
					save `temp6'
						
					rename  `id'	`id0'
					rename  `clustervarname' `clustervarname'0
					rename	`longitude' `longitude0'
					rename  `latitude' `latitude0'	
					rename  `new_point' `new_point'0
					cross using `temp6'
					
					local max_sum_conflicts 1
					while `max_sum_conflicts'>0 {
						cap drop `d'
						tempvar d
						geodist `latitude0' `longitude0'  `latitude' `longitude' , gen(`d')	
						drop if `id0'==`id'
						tempvar conflict
						gen `conflict'=1 if `d' < `between' & `clustervarname'0!=`clustervarname' & `new_point'0==1 & `new_point'==1	& `clustervarname'0!=. & `clustervarname'!=.
						tempvar max_new_point
						bys `id': egen `max_new_point'=max(`new_point')
						tempvar sum_conflicts
						bys `id': egen `sum_conflicts'=sum(`conflict')
						fsum `sum_conflicts'
						local max_sum_conflicts=`r(max)'
						tempvar e
						gen  `e'=1 if `sum_conflicts'==`max_sum_conflicts'
						set seed `seed' 
						tempvar sortorder
						se seed `seed'
						gen double `sortorder' = runiform()		
						tempvar sortorder_id
						bys `id': egen `sortorder_id'=min(`sortorder')
						sort `sum_conflicts' `id' `sortorder'
						tempvar sortorder 
						fegen `sortorder' = group(`sum_conflicts' `sortorder_id') if `e'==1
						fsum `sortorder'
						drop if `sortorder'== `r(max)' & `max_new_point'==1
						keep `id' `clustervarname' `latitude' `longitude' `new_point'
						qui duplicates drop `id', force
						tempfile temp7
						save `temp7'
							
						rename  `id'	`id0'
						rename  `clustervarname' `clustervarname'0
						rename	`longitude' `longitude0'
						rename  `latitude' `latitude0'	
						rename  `new_point' `new_point'0
						cross using `temp7'	
					}	
					keep `id' `clustervarname'  `latitude' `longitude' `new_point'
					qui duplicates drop `id', force
					
					* Finally, drop if distance between two new points larger than `within' (a new point may be at `within' with respect to old points, but not to another newly added point)
					local outside_local 1
					while `outside_local' >0 {
						tempfile temp8
						save `temp8'
							
						rename  `id'	`id0'
						rename  `clustervarname' `clustervarname'0
						rename	`longitude' `longitude0'
						rename  `latitude' `latitude0'	
						rename  `new_point' `new_point'0
						cross using `temp8'	
						
						geodist `latitude0' `longitude0'  `latitude' `longitude' , gen(`d')	
						drop if `id0'==`id'	
						tempvar outside
						gen `outside'=1 if `d' > `within' & `clustervarname'0==`clustervarname' 	& `d'!=.
						fsum `outside' 
						local outside_local=`r(N)'
						tempvar max_outside
						bys `id': egen `max_outside'=max(`outside')
						set seed `seed' 
						tempvar sortorder
						gen double `sortorder' = runiform() if `max_outside'==1
						tempvar max_id_sortorder				
						bys `id': egen double `max_id_sortorder'=max(`sortorder') if `max_outside'==1
						tempvar to_drop
						egen 	`to_drop'=group(`max_outside' `max_id_sortorder') if `max_outside'==1
						fsum     `to_drop'
						if `r(N)'>0 {
							drop if `to_drop'==`r(max)' 
						}	
						keep `id' `clustervarname'  `latitude' `longitude' `new_point'
						qui duplicates drop `id', force
					}		
				}	
			}
		} // closes "if distinct==1 bracket"			
		if `within'==0 {
			gen `clustervarname'=_n
		}
		qui {
			qui drop if `latitude'==. | `longitude'==.
			tempvar sorted_clustervarname
			fegen `sorted_clustervarname'=group(`clustervarname')
			replace `clustervarname'= `sorted_clustervarname'
			drop `sorted_clustervarname'
			label var `clustervarname' "Cluster ID"
			duplicates drop `id', force
			merge 1:1 `latitude' `longitude' using `original0', nogen
			sort `clustervarname'
		}
		if "`no_conflicts'"=="1" {
			noisily di "Note: no conflicts found between eligible clusters"
		}	
	} // ends first qui{
end
	
