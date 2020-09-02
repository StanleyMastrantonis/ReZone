##################################################################
##################################################################
############Northern jarrah forest MINE REZONE MODEL###############
##################################################################
##################################################################

########################PACKAGES##################################
install.packages("raster", dependencies= TRUE)
require(raster)

###################Study Parameters###############################
res = 30  #study resolution
CRS = CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") #study coordinate system
radius = 1000 #set vegetation radius parameter
buffer = 500 #set exponential decay buffer

########################SECTION 1#################################
##################################################################

##################################################################
######################DATA READ IN################################
##########################&#######################################
#######################Functions##################################
##################################################################
setwd("") #set the working directory where the required data is stored

##################################################################
#######################Study Area#################################
study_area = shapefile('Study_area')         #model area  

######################Nest Data###################################
nests = read.csv('nests.csv')[,2:4]
nests = subset(nests, nests$Nest_type =="Nest" 
               | nests$Nest_type =="Possible nest" ) #In this study we are only considering Nests

nests_matrix = as.matrix(nests[,c(2,3)], rownames.force = FALSE)

######################Water Point Data############################
water = read.csv('water.csv')[,2:4]

#assign weights based on water type--------------------------------
water$weight = water_weight = 0
water$weight = water_update = ifelse(water$Type == "Ephemeral",0.75,
                                     ifelse(water$Type == "Permanent",1,water_weight)) #adding weights to water point types

water_matrix = as.matrix(water[,2:4], rownames.force = FALSE)

water_sp =  SpatialPointsDataFrame(coords = data.frame(water$X_UTM,water$Y_UTM),
                                   data = water,proj4string = CRS)
#water points within the study boundary-----------------------------
water_intersect = setNames(na.omit(data.frame(over(water_sp, study_area)[,1],
                                              water$X_UTM,water$Y_UTM)),
                           c('ID','X_UTM','Y_UTM'))[,2:3] # set of water falling within study site only

#####################NDVI Data#################################
NDVI = raster('NDVI.tif')
NDVI[NDVI < 0|is.na(NDVI)] = 0 #only need to consider vegetation

####################Mining Schedule##############################
mine_sched = intersect(buffer(shapefile('mining.shp'),0), study_area)

################################Rehab##########################
rehab = intersect(buffer(shapefile("Rehab"),0), study_area)

################FUNCTIONS#####################################

dist_func = function(x,y){
  distlist = (lapply(1:nrow(x), function(xi){
    cbind(x[xi,1],x[xi,2],xi, (apply(y,1, function(yi){sqrt(sum((x[xi,]-yi)^2))})))
  }))
  distframe = do.call(rbind, distlist)
  colnames(distframe) = c('X','Y','Index','Distance')
  return(distframe)
}

mindist_func = function(x,y){
  (sapply(1:nrow(x), function(xi){ 
    min(apply(y,1, function(yi){sqrt(sum((x[xi,]-yi)^2))}))
  }))
}

whichmindist_func = function(x,y){
  distlist = (lapply(1:nrow(x), function(xi){ 
    which.min(apply(y,1, function(yi){sqrt(sum((x[xi,]-yi)^2))}))
  }))
  distframe = do.call(rbind, distlist)
  return(distframe)
}

whichmindist_func_par = function(x,y){
  clust = makeCluster(detectCores())
  clusterExport(clust, varlist = c("whichmindist_func_par"))
  distlist = (lapply(1:nrow(x), function(xi){ 
    which.min(parApply(clust,y,1, function(yi){sqrt(sum((x[xi,]-yi)^2))}))
  }))
  distframe = do.call(rbind, distlist)
  stopCluster(clust)
  return(distframe)
}

scale_func =  function(x) { (x-min(x))/(max(x)-min(x))}

scale_func_2 =  function(x) { 2*((x-min(x))/(max(x)-min(x)))-1}

sigmoid_func = function(x, coef, d) {
  1 / (1 + exp((x-d)/coef))
}

sigmoid_func_2 = function(x, coef, d, max.add) {
  max.add / (1 + exp((x-d)/coef))
}

weight_mean = function(x){
  y = x[x > 0]
  (sum((rev(sort((x[x > 0])))*(rev(seq(1,length(y),1)))))/(sum(rev(seq(1,length(y),1)))))
}

weight_mean_dim = function(x,length){
  
  if(is.na(sum(x)) | sum(x) <=0 ){
    return(0)
  } else {
    x_sub = x[x > 0]
    x_sort = sort(x_sub, decreasing = TRUE)
    if (length(x_sort) >= length){
      x_lim = x_sort[1:length]
    } else {
      x_lim =  x_sort[1:length(x_sort)]
    }  
    
    y = seq(1,length(x_lim),1)
    ydim = ifelse(y ==1, 1,
                  ifelse(y == 2,0.8,
                         ifelse(y == 3, 0.4,
                                ifelse(y == 4, 0.2,0))))
    ysum =  sum(x_lim*ydim)
    mean = ifelse(yy >1,1,yy)
    return(xy)
  }
}

geo_mean = function(x){prod(x[x > 0])^(1/length(x[x > 0]))}

sum_func = function(x){
  if(is.na(sum(x)) | sum(x) <=0 ){
    return(0)
  } else {
    x_sub = x[x > 0]
    y = sum(x_sub)
    sum = ifelse(y >1,1,y)
    return(sum)
  }
}

find_func = function(seq, vector){
  y = which(abs(seq-vector) == min(abs(seq-vector)))
  y[1]
}

find_func_all = function(seq, vector){
  which(abs(seq-vector) == min(abs(seq-vector)))
}


########################SECTION 2#################################
##################################################################

##################################################################
######################Topologies##################################
##########################&#######################################
##############Grid Topology Data Matching#########################
##################################################################

#retrieve extents of study_area----------------------------------------
bb = res*round(bbox(study_area)/res)

#Define the spatial grid------------------------------------------
gt = GridTopology(cellcentre.offset = bb[,1], 
                  cellsize = c(res, res),
                  cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/res) + 1)
#Create point topologies (might take several minutes)--------------
sp = SpatialPoints(gt, proj4string = CRS)
sp_topo_study = over(sp, study_area)
sp_topo_mining = over(sp, mine_sched)
sp_topo_rehab = over(sp,rehab)


#Point topologies cleanup------------------------------------------ 
topo_coords_study = setNames(cbind(coordinates(sp), sp_topo_study)[!is.na(sp_topo_study[,1]),1:2],c('X_UTM', 'Y_UTM'))
topo_coords_mining = setNames(cbind(coordinates(sp), sp_topo_mining)[!is.na(sp_topo_mining[,1]),1:2],c('X_UTM', 'Y_UTM'))
topo_coords_rehab = setNames(cbind(coordinates(sp), sp_topo_rehab)[!is.na(sp_topo_rehab[,1]),1:2],c('X_UTM', 'Y_UTM'))


#Topolgies to dataframes and matrices ------------------------------
study_grid_matrix = as.matrix(topo_coords_study, rownames.force = FALSE) 
mining_grid = as.matrix(topo_coords_mining, rownames.force = FALSE) 
rehab_grid = cbind(as.matrix(topo_coords_rehab, rownames.force = FALSE), rehab =1) 
###################Data-topology match###############################

#Match location of nests to grid topologies---------------------
nests_grid = whichmindist_func(nests_matrix,study_grid_matrix)
nests_match = study_grid_matrix[nests_grid,]

#remove duplicates--------------------------
nests_nodups = nests_match[-which(duplicated(nests_match[,1:2])), ]

#Match location of water to coordinates of grid------------------

water_grid = whichmindist_func(water_matrix[,1:2],study_grid_matrix)
water_match = study_grid_matrix[water_grid,]

#waterpoints inside of model area------------------

endo_water_grid_ind = whichmindist_func(as.matrix(water_intersect), study_grid_matrix)
endo_water_grid = study_grid_matrix[endo_water_grid_ind,]
endo_water_match = cbind(endo_water_grid,
                         water_matrix[whichmindist_func(endo_water_grid,water_matrix[,1:2]),3])

########################SECTION 3#################################
##################################################################

##################################################################
####################Nests X Food Value############################
##################################################################

food_value_function = function(nests,grid, NDVI){
  #Transform NDVI to coordinates---------------------
  extract_sp =  SpatialPointsDataFrame(coords = grid, data = data.frame(grid),
                                       proj4string = CRS)
  
  extract_NDVI = data.frame(extract(NDVI, extract_sp))
  names(extract_NDVI)[1] = 'NDVI'
  extract_NDVI$X_UTM = grid[,1]
  extract_NDVI$Y_UTM = grid[,2]
  NDVI_coords = as.matrix(extract_NDVI[,2:3])
  
  #calculate distance to radius around nests---------------------
  
  dist_list = lapply(1:nrow(nests), function(xi){
    which(apply(NDVI_coords,1, function(y){
      sqrt(sum((nests[xi,]-y)^2))} <= radius ))
    
  }
  )
  
  #Match nests to NDVI----------------------------
  match_list = setNames(data.frame(do.call(rbind,lapply(1:length(dist_list), function(xi){
    cbind(nests[xi,1],nests[xi,2],xi,extract_NDVI[dist_list[[xi]],1])
  }
  ))),c('X_UTM','Y_UTM','Index','NDVI'))  
  
  #Cumulative NDVI for each nest----------------------------
  NDVI_agg = setNames(aggregate(match_list$NDVI,by=list(match_list$X_UTM,match_list$Y_UTM),sum)
                      ,c('X_UTM','Y_UTM','NDVI'))
  
  
  nests_NDVI = merge(nests,NDVI_agg, by = c('X_UTM',   'Y_UTM'))
  
  #Assign food values to each nest ----------------------------
  food_val = scale_func(seq(min(nests_NDVI$NDVI),max(nests_NDVI$NDVI),1))
  sig_foodval = scale_func(sigmoid_func(food_val, -0.3,-5))
  foodval_df = setNames(data.frame(seq(min(nests_NDVI$NDVI),
                                       max(nests_NDVI$NDVI),1),sig_foodval),c('food_seq','sig_foodval'))
  
  #Match values to data--------------------------------------
  nests_NDVI$foodvalue = foodval_df$sig_foodval[sapply(nests_NDVI$NDVI,find_func_all,foodval_df$food_seq)]
  food_nodups = nests_NDVI[which(!duplicated(nests_NDVI[,1:2])),][,4]
  return(food_nodups)
}

food = food_value_function(nests_nodups, study_grid_matrix,NDVI)

########################SECTION 4#################################
##################################################################

##################################################################
###################Landscape Value Function ######################
##################################################################


maxdist_base = max(mindist_func(nests_nodups,water_matrix[,1:2]))
mindist_base = 0
mindist = mindist_func(nests_nodups,water_matrix[,1:2])

#Value landscapes based on all water points in range and food values---------------------------

landscape_value = function(nests,water,food){
  
  nest_value = sigmoid_func_2(scale_func(seq(mindist_base,maxdist_base,1)), 0.08,0.5,1)
  nest_value_df = setNames(data.frame(seq(mindist_base,maxdist_base,1)
                                      ,nest_value),c('value_seq','nest_value'))
  
  value_list =  setNames(data.frame(do.call(rbind,lapply(1:nrow(nests), function(x){
    cbind(nests[x,1],nests[x,2],x, do.call(rbind, (lapply(1:nrow(water), function(y){
      cbind(sqrt((nests[x,1]-water[y,1])^2+(nests[x,2]-water[y,2])^2), water[y,3]
            
      )}
    ))))
  }))),  c('X_UTM','Y_UTM','hollow','distance', 'weight'))
  
  value_list$inside  = value_list$distance <= max(mindist)
  value_list$value = ind = ifelse(value_list$inside == FALSE,0,
                                  nest_value_df$nest_value[sapply(value_list$distance,
                                                                  find_func,nest_value_df$value_seq)])
  
  value_list$value_prod = value_list$value*value_list$weight
  agg_value = setNames(aggregate(value_list$value_prod,
                                 by=list(value_list$hollow),sum_func), c('hollow','value'))
  landscape_total_value  = sum(agg_value$value*food)
  
  return(landscape_total_value)
  
}

max_val = landscape_value(nests_nodups, endo_water_match, food)


################################################################
#####################Cell value function########################
################################################################
#calculate value of each point across the landscape nests------- 

cellvalue = function(nests, water, food){
  
  maxvalue = landscape_value(nests, water, food)
  #calculate value of pixel for nests--------------------------- 
  cell_value = maxval_new = xxi = yyi = NULL
  for (n in 1:nrow(nests)) {
    nests_ind = which(nests[n,1] == nests[,1]
                      & nests[n,2] == nests[,2])
    nests_mat = nests[-nests_ind,]
    food_ind = food[-nests_ind]
    maxvalue_new =landscape_value(nests_mat, water, food_ind)
    cell_value[n] = maxvalue - maxvalue_new
    xxi[n] = nests[n,1]
    yyi[n] = nests[n,2]
  }
  
  nest_value = data.frame(xxi,yyi,cell_value)
  cell_value = xxi = yyi = NULL
  
  for (m in 1:nrow(water)){
    water_ind = which(water[m,1] == water[,1]
                      & water[m,2] == water[,2])
    water_matrix = water[-water_ind,]
    maxvalue_new = landscape_value(nests, water_matrix, food)
    cell_value[m] = ifelse((maxvalue - maxvalue_new) > 1,1, maxvalue - maxvalue_new)
    xxi[m] = water[m,1]
    yyi[m] = water[m,2]
  }
  
  water_value = data.frame(xxi,yyi,cell_value)
  study_value = rbind(nest_value, water_value)
  return(study_value)
  
}


cell_value = cellvalue(nests_nodups, endo_water_match, food)


################################################################
####################Value Buffer Function#######################
################################################################

#This function will return a raster with values using exponential decay------- 

raster_value_func = function(cellvalue, grid){
  
  buff_seq = seq(1,buffer,1)
  y_range = -1*(buff_seq^0.2) + buffer
  y_scale = scale_func(y_range)
  buffer_df = data.frame(buff_seq,y_scale)
  cellvalue_mat = as.matrix(cellvalue[,1:2])
  
  distbuffer = lapply(1:nrow(cellvalue_mat), function(xi)
    which(apply(grid,1, function(y)
      sqrt(sum((cellvalue_mat[xi,]-y)^2)) <= buffer ))
  )
  
  match_buffer = setNames(data.frame(do.call(rbind,lapply(1:length(distbuffer), function(xi){
    cbind(grid[distbuffer[[xi]],1],grid[distbuffer[[xi]],2],xi,
          sapply(1:length(distbuffer[[xi]]), function(yi){
            sqrt(sum((cellvalue_mat[xi,1]-grid[distbuffer[[xi]][yi],1])^2)+
                   (cellvalue_mat[xi,2]-grid[distbuffer[[xi]][yi],2])^2)}
          ))
  }))),c('X_UTM','Y_UTM','index','distance')) 
  
  match_buffer$value= buffer_df$y_scale[sapply(match_buffer$distance,
                                               find_func_all,buffer_df$buff_seq)]
  
  agg_buff = setNames(aggregate(match_buffer$value,
                                by=list(match_buffer$X_UTM, match_buffer$Y_UTM),sum_func),c('X_UTM','Y_UTM','value') )
  
  merge_buffer = merge(grid, agg_buff, by = c('X_UTM','Y_UTM'), all.x = TRUE)
  merge_buffer$value = ifelse(is.na(merge_buffer$value),0,merge_buffer$value)
  raster_value = rasterFromXYZ(merge_buffer, crs = CRS)
  
  return(raster_value)
  
}  


ras_value = raster_value_func(cell_value, study_grid_matrix)
grid_value = data.frame(rasterToPoints(ras_value))
colnames(grid_value) = c('X_UTM','Y_UTM','value')

#######you can save this raster:
#writeRaster(ras_value, 'value_ras', 'GTiff', prj = TRUE)

########################SECTION 5#################################
##################################################################

##################################################################
################Pareto optimal Frontier###########################
##################################################################

#######This algorithm will generate the pareto effeciency frontier and save the landscape states####
mining_grid_df = data.frame(mining_grid)
mining_grid_df$mine = 1
mine_merge = merge(grid_value,mining_grid_df, by = c('X_UTM','Y_UTM'), all.x = TRUE)
mine_merge$mine = mine_check = ifelse(is.na(mine_merge$mine) == TRUE, 0,mine_merge$mine )

mine_seq_loop = seq(min(mine_merge$value),max(mine_merge$value), 0.01)
nomine_baseline = sum(mine_merge$value)
mine_extent = mine_amount = value_mine = NULL
landstate_df = setNames(data.frame(matrix(ncol = 4, nrow = 0)), names(mine_merge))
subset_mine_0 = subset(mine_merge, mine_merge$mine == 0)
mine_0_sum = sum(subset_mine_0$mine)
value_mine[1] = mine_extent[1] = value_calc =  nomine_baseline
mine_amount[1] = mine_0_sum
range_new = value_range = range_count = 0  
counter_1 = 2
counter_2 = 1

sub_choice = mine_amount_ls = list()

for (range_count in mine_seq_loop){
  value_calc_new = value_calc
  
  if (range_count == 0){
    subset_landscape = subset(mine_merge, mine_merge$value == 0 
                              & mine_merge$mine == 1 )
  } else {
    subset_landscape = subset(mine_merge, mine_merge$value <= range_count &
                                mine_merge$value > range_new &
                                mine_merge$mine == 1 )
  }
  
  landstate_df = rbind(landstate_df,subset_landscape)
  mine_sum = sum(subset_landscape$mine)
  value_sum = sum(subset_landscape$value)
  range_new = range_count
  value_range = value_sum
  value_calc = value_calc_new - value_range
  value_mine[counter_1] = value_calc
  mine_extent = mine_amount[counter_2]
  mine_amount[counter_1] = mine_sum+mine_extent
  frontier_value_ind = mine_amount
  sub_choice[[counter_2]] = landstate_df
  mine_amount_ls[[counter_2]] = mine_amount
  counter_1= counter_1+1
  counter_2= counter_2+1
  
}


df_frontier = data.frame(value_mine,mine_amount)
df_frontier$percentage = df_frontier$value_mine/nomine_baseline*100
df_frontier$area = df_frontier$mine_amount*900/1000000


options(scipen = 999)
plot(df_frontier$area,df_frontier$percentage, main = 'Pareto Efficiency Frontier', ylab ='Total Landscape Value (%)',
     xlab = expression(paste("Area of Mining"," (km"^"2",')')), type = 'o', col = 'red', lwd = 4, cex = 0)

#######IMPORTANT#######################################################
#######click on the fontier plot to choose a land use allocation 'finish' in top right of plot
choice = identify(df_frontier$area, 
                  df_frontier$percentage, pos = FALSE, tolerance = 5 )
##########################################################################
v = df_frontier[choice,4]
abline(v =  v, lty = 5, lwd = 2)
legend("topright", legend=c("Frontier", "User Decision"),
       col=c("red", "black"), lty=c(1,5),lwd = 2, cex=1)

choice_ls = sub_choice[[choice-1]][,c(1,2,4)]
choice_ras = rasterFromXYZ(choice_ls, crs = CRS)
plot(choice_ras, legend = FALSE)

#######you can save this raster:
#writeRaster(choice_ras, 'choice_ras', 'GTiff', prj = TRUE)
########################SECTION 6#################################
##################################################################

##################################################################
################Heuristic Algorithm###############################
##################################################################

###################This might take some hours if you let it run fully################
max_val_ini = landscape_value(nests_nodups, endo_water_match, food)
max_val_new = max_val_ini
orthogonal_1 = res+1  
diagonal_1 = sqrt(res^2+res^2)+1
orthogonal_2 = res*2+1
diagonal_2 = sqrt((res*2)^2+(res*2)^2)+1
orthogonal_3 = res*3+1
diagonal_3 = sqrt((res*3)^2+(res*3)^2)+1

iteration = 1
counter = 1
recurse = 'FALSE'

confirmed_points = matrix(ncol =3, nrow = 0) 
colnames(confirmed_points) = c(' X_UTM',   'Y_UTM', 'count')
failed_point = matrix(ncol =2, nrow = 0) 
break_count =1
num_water_points = 10 
plot(NDVI)
points(nests_nodups, col = 'red', pch = 20)
points(water_matrix, col = 'light blue', pch = 20, cex = 2)
legend('topleft', legend=c("Nests", "Water"),
       col=c("red", "lightblue"), pch = 20, cex=1)
count = 1

while (count <= num_water_points){
  
  break_count = break_count +1
  if (break_count > num_water_points+10){
    break
  }
  
  if (count == 1){
    newmax_prev = max_val_ini
    iteration = 1
    newmax_vec = vector('numeric')
    water_mat = water_matrix
    newpoint_matrix = matrix(ncol =2, nrow =0)
    colnames(newpoint_matrix) = c('X_UTM','Y_UTM')
    recurse == 'FALSE'
    old_vec = vector('numeric')
    old_max = max_val_ini
    
  } else {
    
    water_mat = water_mat
    newmax_prev = newmax_prev
    iteration = 1
    newpoint_matrix = newpoint_matrix
    recurse == 'FALSE'
    old_vec = vector('numeric')
    old_max = newmax_prev
    
  }
  
  while(iteration <= 3){
    
    if (iteration == 1 & recurse == "FALSE") {
      diag = diagonal_1
      ortho = orthogonal_1
      min_diag = 0
      min_ortho = 0
      rand_point = as.matrix(t(study_grid_matrix[sample(nrow(study_grid_matrix), 1), ]),
                             rownames.force = FALSE)
      land_matrix = study_grid_matrix
      newmax_prev = newmax_prev
      
    } else if (iteration == 1 & recurse == 'TRUE'){
      
      rand_point = new_point
      diag = diagonal_1
      ortho = orthogonal_1
      min_diag = 0
      min_ortho = 0
      land_matrix = study_grid_matrix[-old_vec,1:2]
      row.names(land_matrix) = NULL
      newmax_prev = newmax_prev
      
    } else if (iteration == 2){
      
      rand_point = new_point
      diag = diagonal_2
      ortho = orthogonal_2
      min_diag = diagonal_2
      min_ortho = orthogonal_2
      land_matrix = study_grid_matrix[-old_vec,1:2]
      row.names(land_matrix) = NULL
      newmax_prev = newmax_prev
      
    } else if (iteration == 3){
      
      rand_point = new_point
      diag = diagonal_3
      ortho = orthogonal_3
      min_diag = diagonal_2
      min_ortho = orthogonal_2
      land_matrix = study_grid_matrix[-old_vec,1:2]
      row.names(land_matrix) = NULL
      newmax_prev = newmax_prev
      
    }
    
    points(rand_point, col = 'black', pch = 17, cex = 1.5)
    old_point = rand_point
    which_old = as.matrix(t(study_grid_matrix[which(old_point[1,1] == study_grid_matrix[,1] &
                                                      old_point[1,2] == study_grid_matrix[,2]),]))
    old_vec = rbind(old_vec,which_old)
    row.names(old_vec) = NULL
    distmat_wp = apply(land_matrix,1, function(x){sqrt(sum((x-rand_point)^2))})
    land_matrix = cbind(land_matrix, distmat_wp)
    colnames(land_matrix)[3] = 'dist'
    land_matrix = land_matrix[order(land_matrix[,3]),]
    land_matrix[,3] = as.numeric(format(round(land_matrix[,3], 5), nsmall = 5))
    which_mat = which(land_matrix[,3] > min_ortho & land_matrix[,3] <=  diag)
    land_matrix = land_matrix[which_mat,]
    val_wp = NULL
    new_wp = NULL
    new_wp_bind = NULL
    
    for ( i in 1:nrow(land_matrix)){
      
      new_wp = cbind(as.matrix(t(land_matrix[i,1:2])),weight = 1.00)
      new_wp_bind = rbind(water_mat,new_wp)
      row.names(new_wp_bind) = NULL
      val_wp[i] = landscape_value(nests_nodups,new_wp_bind, food_nodups)
      
    }
    
    val_wp = val_wp[!is.na(val_wp)]
    wp_ind = which(val_wp == max(val_wp))
    wp_ind = wp_ind[1]
    newmax = val_wp[wp_ind]
    print(newmax)
    
    if (newmax > newmax_prev){
      iteration = 1
      recurse = 'TRUE'
      new_point = as.matrix(t(land_matrix[wp_ind,1:2]))
      newmax_vec = append(newmax_vec, newmax)
      newpoint_matrix = rbind(newpoint_matrix , new_point)
      newmax_prev = max(newmax_vec)
      points(new_point, col = 'black', pch = 17, cex = 1.5)
      
      
    } else {
      
      iteration = iteration + 1
      recurse = 'TRUE'
      new_point = old_point
      newmax_prev = newmax_prev
      
    }
    
    if (iteration > 3){
      
      recurse = 'FALSE'
      iteration = 1
      break
      
    }
  }
  
  if (count == 1){
    
    if (newmax > newmax_prev){
      
      new_point_count = new_point
      new_point_count = cbind(new_point, count)
      confirmed_points = rbind(confirmed_points, new_point_count)
      count = count + 1
      
    } else {
      
      if (newmax > old_max){
        
        new_point_count = new_point
        new_point_count = cbind(new_point, count)
        confirmed_points = rbind(confirmed_points, new_point_count)
        count = count + 1
        
      } else {
        
        failed_point = rbind(failed_point, new_point)
        count = count 
        
      }
    }
    
  } else {
    
    if (newmax > old_max){
      
      new_point_count = new_point
      new_point_count = cbind(new_point, count)
      confirmed_points = rbind(confirmed_points, new_point_count)
      count = count + 1
      
    } else {
      
      failed_point = rbind(failed_point, new_point)
      
    }
  }
  
  max_ind = which(newmax_vec == max(newmax_vec))
  newmax_prev = newmax_vec[max_ind]
  new_water_point = cbind(as.matrix(t(newpoint_matrix[max_ind,])),weight = 1)
  water_mat = rbind(water_mat,new_water_point)
  row.names(water_mat) = NULL
  print(paste('Count = ', count))
  
}



water_mat #####this is the result of 10 new water points


#################################################################
#################################Zoning##########################


cons_value = grid_value
cons_grid = merge(data.frame(study_grid_matrix),cons_value, all.x = TRUE )
cons_grid$value = val_zone = ifelse(is.na(cons_grid$value),0,cons_grid$value)

rehabilitation_grid = merge(cons_grid,rehab_grid, all.x = TRUE)
rehabilitation_grid$rehab = rehab_zone = ifelse(is.na(rehabilitation_grid$rehab)
                                                ,0,rehabilitation_grid$rehab)

mine_grid = merge(rehabilitation_grid,choice_ls, all.x = TRUE)
mine_grid$mine = mine_zone = ifelse(is.na(mine_grid$mine),0,mine_grid$mine)

zone_grid = mine_grid

zone_grid$class = classes = as.factor(ifelse(zone_grid$value > 0.10, 'Conservation',
                                             ifelse(zone_grid$rehab == 1, 'Rehabilitation',
                                                    ifelse(zone_grid$mine == 1, 'Mining','NA'))))

zone_grid$int = class_int = ifelse(classes == 'Conservation', 1,
                                   ifelse(classes ==  'Rehabilitation',2,
                                          ifelse(classes == 'Mining',3,4)))

zone_df = zone_grid[,c(1:2,7)]



zone_ras = rasterFromXYZ(zone_df, crs = CRS)
zone_att = ratify(zone_ras)
rat_lv = levels(zone_att)[[1]]
rat_lv$class_name = levels(zone_grid$class)
levels(zone_att) = rat_lv
zoning = deratify(zone_att, 'class_name')
plot(zone_att)
#######you can save this raster:
#writeRaster(zoning, 'zone_ras', 'GTiff', prj = TRUE , RAT=TRUE )