cppFunction('NumericMatrix crossdist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();
  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }
  NumericMatrix out(nrow1, nrow2);
  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1,r2) = sqrt(total);
    }
  }
  return out;
}')

cppFunction('NumericVector pairdist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int ncol = m1.ncol();
  NumericVector out(nrow1);
  for (int r1 =0 ; r1 < nrow1; r1++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r1, c12), 2);
      }
      out(r1) = sqrt(total);
    
  }
  return out;
}')



mindist_func = function(x,y){
  d_list = list()
  cd = crossdist(as.matrix(x),as.matrix(y) )
  for(i in 1:ncol(cd)) {           
    d_list[[i]] = cd[, i]
  }
  mind = unlist(lapply(d_list, min))
  return(mind)
}

maxdist_func = function(x,y){
  d_list = list()
  cd = crossdist(as.matrix(x),as.matrix(y) )
  for(i in 1:ncol(cd)) {           
    d_list[[i]] = cd[, i]
  }
  mind = unlist(lapply(d_list, max))
  return(mind)
}


whichmindist_func = function(x,y){
  d_list = list()
  cd = crossdist(as.matrix(x),as.matrix(y) )
  for(i in 1:ncol(cd)) {           
    d_list[[i]] = cd[, i]
  }
  # print(d_list)
  which_min = unlist(lapply(d_list, which.min))
  return(which_min)
}

whichdist_func = function(x,y,dist){
  d_list = list()
  cd = crossdist(as.matrix(x),as.matrix(y) )
  for(i in 1:ncol(cd)) {           
    d_list[[i]] = cd[, i]
  }
  which = unlist(lapply(d_list, function(x){which(x <= dist )} ))
  return(which)
}

whichdist_func_list = function(x,y,dist){
  d_list = list()
  cd = crossdist(as.matrix(x),as.matrix(y) )
  for(i in 1:ncol(cd)) {           
    d_list[[i]] = cd[, i]
  }
  which = lapply(d_list, function(x){which(x <= dist )} )
  return(which)
}




scale_func =  function(x) { (x-min(x))/(max(x)-min(x))}

scale_func_range =  function(x, a, b) { (b-a)*((x-min(x))/(max(x)-min(x)))+a}

scale_func_ras = function(x){
  (x-(cellStats(x, stat='min')))/((cellStats(x, stat='max'))-(cellStats(x, stat='min')))
}


scale_func_range_ras = function(x, a ,b){
  (b-a)*(x-(cellStats(x, stat='min')))/((cellStats(x, stat='max'))-(cellStats(x, stat='min')))+a
}


sigmoid_func = function(x, coef, d) {
  1 / (1 + exp((x-d)/coef))
}

sigmoid_func_2 = function(x, coef, d, max.add) {
  max.add / (1 + exp((x-d)/coef))
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

mean_func = function(x){
  if(is.na(sum(x)) | sum(x) <=0 ){
    return(0)
  } else {
    x_sub = x[x > 0]
    y = mean(x_sub)
    return(y)
  }
  
}


find_func = function(seq, vector){
  y = which(abs(seq-vector) == min(abs(seq-vector)))
  y[1]
}

find_func_all = function(seq, vector){
  which(abs(seq-vector) == min(abs(seq-vector)))
}



idw_func_k = function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    x_len = ifelse(length(weight) > max_kern,max_kern, (length(weight)))
    xw = weight/distance
    xw_x = sum(xw)
    pd = sum(1/distance)
    dw = (xw_x/pd)*(x_len/max_kern)
    return(dw)
  }
  
}


d_sum_func = function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    xw = 1-(weight/(distance+weight))
    xw_x = 1-(prod(xw))
    return(xw_x)
  }
  
}

d_sum_func_2 = function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    xw = 1-(1-weight)
    xw_x = 1-(prod(xw))
    return(xw_x)
  }
  
}

idw_func = function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    
    xw = weight/distance
    # ind = which(is.infinite(xw))
    xw_x = sum(xw)
    pd = sum(1/distance)
    dw = (xw_x/pd)
    return(dw)
  }
  
  
}


area_func = function(patch){
  # if(nrow(patch) < 1){
  #   return(0)
  # } else {
    
    a_p = sum(patch[patch>0])
    print(a_p)
    area = length(patch) 
    # ind = which(is.infinite(xw))
    prop = a_p/area
    return(prop)
  # }
  
  
}


get_topologies = function(area, res){
    bb = res*round(bbox(study_area)/res)
    gt = GridTopology(cellcentre.offset = bb[,1], 
                      cellsize = c(res, res),
                      cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/res) + 1)
    
    sp = SpatialPoints(gt)
    sp_topo_study = over(sp, study_area)
    topo_coords_study = data.frame(cbind(coordinates(sp), sp_topo_study)[!is.na(sp_topo_study),1:2])
    
    grid = as.matrix(topo_coords_study, rownames.force = FALSE)
    return(grid)
}

sdf = function(min,max){
  maxdist_base = max
  mindist_base = min
  value = sigmoid_func_2(scale_func(seq(mindist_base,maxdist_base,0.1)), 0.08,0.5,1)
  value_df = data.frame(seq(mindist_base,maxdist_base,0.1),value)
  colnames(value_df) = c('value_seq','value')
  return(value_df)
}



max_kern_func = function(resource1, resource2){
  kern_d_list = data.frame(matrix(crossdist(resource1, resource2), ncol =1))
  kern_d_list$index = rep(seq(1,nrow(resource1),1), times = nrow(resource2))
  colnames(kern_d_list) =  c('distance','index' )
  kern_d_list$inside  = kern_d_list$distance <= forage
  kern_list = kern_d_list[kern_d_list$inside == TRUE,]
  max_kern = ceiling(median(table(kern_list$index)))
  return(max_kern)
}



res_value_function = function(resource1,grid, ac){
  
  ac_mat = rasterToPoints(ac)
  colnames(ac_mat) = c('x','y','val')
  ac_coords = as.matrix(ac_mat[,1:2])
  dist_list = whichdist_func_list(ac_coords, resource1, forage)
  match_list = data.frame(do.call(rbind,lapply(1:length(dist_list), function(xi){
    cbind(resource1[xi,1],resource1[xi,2],xi,ac_mat[dist_list[[xi]],3])
  }
  )))
  
  match_list$ac_ind =unlist(dist_list)
  colnames(match_list) = c('x','y','index','ac','ac_ind')
  match_list$ac_x = ac_coords[match_list$ac_ind,1]
  match_list$ac_y = ac_coords[match_list$ac_ind,2]
  match_list$dist = pairdist(as.matrix(match_list[,1:2]),as.matrix(match_list[,6:7] ))
  ind_loop = unique(match_list$index)
  idw_val = NULL
  x = NULL
  y = NULL
  
  for (i in 1:length(ind_loop)){
    sub = match_list[match_list$index == ind_loop[i] & match_list$dist > 0,]
    idw_val[i] = ifelse(is.null(dim(sub)),0,idw_func(sub$ac,sub$dist))
    x[i] = sub[1,1]
    y[i] = sub[1,2]
  }
  
  # df_agg = data.frame(x,y,scale_func(idw_val))
  df_agg = data.frame(x,y,idw_val)
  colnames(df_agg) =  c('x','y','ac')
  # return(df_agg)
  return(df_agg$ac)
  
}


res_value_function_area = function(resource1,grid, ac){
  
  ac_mat = rasterToPoints(ac)
  colnames(ac_mat) = c('x','y','val')
  ac_coords = as.matrix(ac_mat[,1:2])
  dist_list = whichdist_func_list(ac_coords, resource1, forage)
  match_list = data.frame(do.call(rbind,lapply(1:length(dist_list), function(xi){
    cbind(resource1[xi,1],resource1[xi,2],xi,ac_mat[dist_list[[xi]],3])
  }
  )))
  
  match_list$ac_ind =unlist(dist_list)
  colnames(match_list) = c('x','y','index','ac','ac_ind')
  match_list$ac_x = ac_coords[match_list$ac_ind,1]
  match_list$ac_y = ac_coords[match_list$ac_ind,2]
  match_list$dist = pairdist(as.matrix(match_list[,1:2]),as.matrix(match_list[,6:7] ))
  ind_loop = unique(match_list$index)
  area_val = NULL
  x = NULL
  y = NULL
  
  for (i in 1:length(ind_loop)){
    sub = match_list[match_list$index == ind_loop[i] & match_list$dist > 0,]
    area_val[i] = ifelse(is.null(dim(sub)),0,area_func(sub$ac))
    x[i] = sub[1,1]
    y[i] = sub[1,2]
  }
  
  # df_agg = data.frame(x,y,scale_func(area_val))
  df_agg = data.frame(x,y,area_val)
  colnames(df_agg) =  c('x','y','ac')
  # return(df_agg)
  return(df_agg$ac)
  
}

landscape_value = function(resource1,resource2,ac_val){
  
  value_list = data.frame(matrix(crossdist(resource1,resource2), ncol =1))
  value_list$index = rep(seq(1,nrow(resource1),1), times = nrow(resource2))
  colnames(value_list) =  c('distance','index' )
  value_list$inside  = value_list$distance <= forage
  value_list$value  = value_df$value[sapply(value_list$distance,
                                            find_func,value_df$value_seq)]
  ind_loop = unique(value_list$index)
  idw_val = NULL
  
  for (i in 1:length(ind_loop)){
    sub = value_list[value_list$index == ind_loop[i],]
    idw_val[i] = ifelse(all(sub$inside == FALSE) ,0,idw_func(sub[sub$inside == TRUE, 4],sub[sub$inside == TRUE, 1]))
  }
  
  agg_value = data.frame(resource1[,1],resource1[,2],idw_val)
  colnames(agg_value) = c('x','y','value')
  landscape_total_value  = sum(agg_value$value*ac_val)
  return(landscape_total_value)
  
}

landscape_value_sum = function(resource1,resource2,ac_val){
  
  value_list = data.frame(matrix(crossdist(resource1,resource2), ncol =1))
  value_list$index = rep(seq(1,nrow(resource1),1), times = nrow(resource2))
  colnames(value_list) =  c('distance','index' )
  value_list$inside  = value_list$distance <= forage
  value_list$value  = value_df$value[sapply(value_list$distance,
                                            find_func,value_df$value_seq)]
  ind_loop = unique(value_list$index)
  sum_val = NULL
  
  for (i in 1:length(ind_loop)){
    sub = value_list[value_list$index == ind_loop[i],]
    sum_val[i] = ifelse(all(sub$inside == FALSE) ,0,d_sum_func_2(sub[sub$inside == TRUE, 4],sub[sub$inside == TRUE, 1]))
  }
  
  agg_value = data.frame(resource1[,1],resource1[,2],sum_val)
  colnames(agg_value) = c('x','y','value')
  landscape_total_value  = sum((agg_value$value*ac_val))
  # return(landscape_total_value)
  return(landscape_total_value)
  
}





cellvalue = function(resource1, resource2, ac_val){
  maxvalue = landscape_value_sum(resource1, resource2, ac_val)
  cell_value = maxval_new = xxi = yyi = NULL
  for (n in 1:nrow(resource1)) {
    res1_new = resource1[-n,]
    val_ind = ac_val[-n]
    maxvalue_new =landscape_value_sum(res1_new , resource2, val_ind)
    cell_value[n] = maxvalue - maxvalue_new
    xxi[n] = resource1[n,1]
    yyi[n] = resource1[n,2]
    
  }
  res1_value = data.frame(xxi,yyi,(cell_value))
  colnames(res1_value) = c('x','y','val')
  cell_value = xxi = yyi = NULL
  
  for (m in 1:nrow(resource2)){
    
    res2_new = resource2[-m,]
    maxvalue_new = landscape_value_sum(resource1, res2_new, ac_val)
    # cell_value[m] = abs(ifelse((maxvalue - maxvalue_new) > 1,1, maxvalue - maxvalue_new))
    cell_value[m] = ifelse((maxvalue - maxvalue_new) > 1 ,1, maxvalue - maxvalue_new)
    xxi[m] = resource2[m,1]
    yyi[m] = resource2[m,2]
    
  }
  
  res2_value = data.frame(xxi,yyi,(cell_value))
  colnames(res2_value) = c('x','y','val')
  study_value = rbind(res1_value, res2_value)
  return(study_value)
  
}


raster_value_func = function(cellvalue, grid){
  
  buff_seq = seq(1,buffer,1)
  y_range = -1*(buff_seq^0.2) + buffer
  y_scale = scale_func(y_range)
  buffer_df = data.frame(buff_seq,y_scale)
  cellvalue_mat = as.matrix(cellvalue[,1:2])

  
  distbuffer = data.frame(matrix(crossdist(grid, cellvalue_mat), ncol =1))
  distbuffer$index = rep(1:nrow(cellvalue_mat), each = nrow(grid))
  distbuffer$value = cellvalue[distbuffer$index,3]
  distbuffer$grid_index = rep(1:nrow(grid), times = nrow(cellvalue))
  distbuffer$x = grid[distbuffer$grid_index,1]
  distbuffer$y = grid[distbuffer$grid_index,2]
  colnames(distbuffer) = c('dist','index','value','grid_index','x','y')
  distbuffer$inside = distbuffer$dist <= buffer
  distbuffer$val_in = ifelse(distbuffer$inside == FALSE,0,distbuffer$value)
  distbuffer = distbuffer[distbuffer$inside == TRUE,]
  distbuffer$weight = buffer_df$y_scale[sapply(distbuffer$dist,
                                               find_func_all,buffer_df$buff_seq)]
  
  distbuffer$w_val = distbuffer$val_in*distbuffer$weight
  agg_buff = aggregate(distbuffer$w_val, list(distbuffer$x, distbuffer$y), sum_func)
  
  
  ras = rasterFromXYZ(agg_buff)
  ras[is.na(ras)] = 0 
  return(ras)
}  


passage_func = function(res1, res2,value, rasval,  perm){
  p_1 = value[1:nrow(res1),]
  p_2 = value[nrow(res1): nrow(value),]
  
  p_1_samp = p_1[sample(1:nrow(p_1), size = 20, prob  = p_1[,3]),]
  # p_2_samp = p_2[sample(1:nrow(p_2), size = 10, prob  = p_2$cell_value),]
  p_2_samp = p_2
  tran = transition(perm, function(x) 1/mean(x), direction =16, symm =FALSE)
  tran_g = geoCorrection(tran)
  
  lraster = c()
  #nrow(p_1)
  for(i in 1:nrow(p_1)){
    opoint = as.matrix(p_1[sample(1:nrow(p_1), size = 1, prob  = p_1[,3]),][,1:2], ncol = 2)
    gpoint = as.matrix(p_2[sample(1:nrow(p_2), size = 1, prob  = p_2[,3]),][,1:2], ncol = 2)
    passraster = passage(tran_g, opoint, gpoint, 0.5)
    lraster = append(lraster, passraster)
    print(i)
  }
  
  sum_rast = Reduce("+",lraster)
  p_v = value
  p = SpatialPointsDataFrame(coords = p_v[,1:2],
                             data = p_v)
  # p$cell_value_2 = 1
  sum_rast[cellFromXY(sum_rast, p@coords)] = 0
  #rasval[cellFromXY(rasval, p@coords)]  
  
  pass_ras = scale_func_ras(sum_rast)
  
  ras_list = stack(pass_ras, rasval)
  ras_sum = calc(ras_list , sum)
  ras_sum[ras_sum > 1] = 1
  # ras_sum[cellFromXY(ras_sum, p@coords)] = rasval[cellFromXY(rasval, p@coords)]
  fix_points = as.matrix(cbind(value[,1:2], ras_sum[cellFromXY(ras_sum, value[,1:2])]))
  colnames(fix_points) = c('x','y','value')
  ras_points = as.matrix(rasterToPoints(ras_sum))
  
  fix_list = data.frame(matrix(crossdist(ras_points[,1:2],fix_points[,1:2]), ncol =1))
  fix_list$index = rep(seq(1,nrow(fix_points),1), each = nrow(ras_points))
  fix_list$val = rep(ras_points[,3], times = nrow(fix_points))
  colnames(fix_list) =  c('distance','index', 'val')
  fix_list$inside  = fix_list$distance <= 2
  fix_list_sub = fix_list[fix_list$inside == TRUE & fix_list$distance > 0,]
  fix_agg = aggregate(fix_list_sub$val, list(fix_list_sub$index), mean)
  ras_sum[cellFromXY(ras_sum, value[,1:2])] = fix_agg$x
  
  return(ras_sum)
  
}


pareto_func = function(ras,conf){
ras = ras
conf[is.na(conf)] = 0
ras[is.na(ras)] = 0

stack(conf,ras)

conf_grid = rasterToPoints(conf)
val_grid = rasterToPoints(ras)
colnames(val_grid) = c('x','y','pass')

merge_land = merge(val_grid,conf_grid, by = c('x','y'), all.x = TRUE)
seq_loop = seq(min(merge_land$pass),max(merge_land$pass), 0.01)


noconf_baseline = sum(merge_land$pass)
conf_baseline = sum(merge_land$conflict)
conf_extent = conf_amount = value_conf = NULL
landstate_df = setNames(data.frame(matrix(ncol = 4, nrow = 0)), names(merge_land))
subset_conf_0 = subset(merge_land, merge_land$conflict == 0)
conf_0_sum = sum(subset_conf_0$conflict)
value_conf[1] = conf_extent[1] = value_calc =  noconf_baseline
conf_amount[1] = conf_0_sum
range_new = value_range = range_count = 0  
counter_1 = 2
counter_2 = 1

sub_choice = conf_amount_ls = list()

for (range_count in seq_loop){
  value_calc_new = value_calc
  
  if (range_count == 0){
    subset_landscape = subset(merge_land, merge_land$pass == 0 
                              & merge_land$conflict == 1 )
  } else {
    subset_landscape = subset(merge_land, merge_land$pass <= range_count &
                                merge_land$pass > range_new &
                                merge_land$conflict == 1 )
  }
  
  landstate_df = rbind(landstate_df,subset_landscape)
  conf_sum = sum(subset_landscape$conflict)
  value_sum = sum(subset_landscape$pass)
  range_new = range_count
  value_range = value_sum
  value_calc = value_calc_new - value_range
  value_conf[counter_1] = value_calc
  conf_extent = conf_amount[counter_2]
  conf_amount[counter_1] = conf_sum+conf_extent
  frontier_value_ind = conf_amount
  sub_choice[[counter_2]] = landstate_df
  conf_amount_ls[[counter_2]] = conf_amount
  counter_1= counter_1+1
  counter_2= counter_2+1
  
}


df_frontier = data.frame(value_conf,conf_amount)
df_frontier$percentage = df_frontier$value_conf/noconf_baseline*100
df_frontier$area = df_frontier$conf_amount/conf_baseline*100


options(scipen = 999)
# par(mfrow = c(4,2), oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 3) + 0.1)
plot(df_frontier$area,df_frontier$percentage, ylab ='Total Habitat Value (%)',
     type = 'o',xlab = "", col = 'red', lwd = 4, cex = 0)
title(xlab = "Area Cleared (%)")
title('Pareto Efficiency Frontier', line = 0.2)

#######IMPORTANT#######################################################
#######click on the fontier plot to choose a land use allocation 'finish' in top right of plot
choice = identify(df_frontier$area,
                  df_frontier$percentage, pos = FALSE, tolerance = 5,
                  atpen = FALSE, labels = '')


##########################################################################
v = df_frontier[choice,4]
abline(v =  v, lty = 5, lwd = 2)
legend("topright", legend=c("Pareto Frontier", "Decision"),
       col=c("red", "black"), lty=c(1,5),lwd = 2, cex=1)

choice_ls = sub_choice[[choice-1]][,c(1,2,4)]
choice_ras = rasterFromXYZ(choice_ls)


# up_ras = ras
# # plot(up_ras)
# 
# up_ras[cellFromXY(up_ras, choice_ls[,1:2])] = 0
# # plot(up_ras)
# 
# plot(up_ras, col = cols_out, legend = TRUE,
#      asp =0, axes = T,xlab = 'x' , interpolate = TRUE,
#      legend.width=2, legend.shrink= 1, legend.args = list(text = '', side = 4, 
#                                                           font = 2, line = 2.5, cex = 1))
# #Habitat Value
# # title('Optimised Conservation Outcome', line = 0.2)
# title(ylab = 'y')
# 
# plot(choice_ras, legend = FALSE, col= 'red', add = TRUE)
# legend(1,300, legend = 'Clearing Extent', col = 'red', pch = 15, cex = 1.3)

return(choice_ras)
}

