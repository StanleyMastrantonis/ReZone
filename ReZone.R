# libraries -------------------
require(raster)
require(Rcpp)
require(gdistance)


# Get the current directory path -------------------
current_dir = getwd()

# Set the working directory to the current directory -------------------
setwd(current_dir)

#functions-------------------
source('ReZone_Functions.R')#source the Rezone functions, please contact author for the heuristic function


#study parameters------------
res = 1  #study resolution
# CRS = CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") #study coordinate system
forage = 75 #set vegetation radius parameter
buffer = 25 #set exponential decay buffer
breaks = seq(0, 1, by = 0.01)
rbPal = colorRampPalette(c('red','blue'))
cols_out = colorRampPalette(c('white',"pink", "yellow", "lime green"))(length(breaks) - 1)

#data read-in----------------
data_dir = file.path(getwd(), "data")
setwd(data_dir)#set directory of data

res_1_mat = as.matrix(read.csv('p_res_1.csv'), rownames.force = FALSE)
res_2_mat = as.matrix(read.csv('p_res_2.csv'), rownames.force = FALSE)
ac = raster('sp_1_ac.tif')
study_area = as(extent(ac), 'SpatialPolygons')
# grad = raster('sp_1_grad.tif')
perm = 1-ac
conf = raster('')#load in conflicting spatial layer

grid = get_topologies(study_area, res)
res_1_match = grid[whichmindist_func(grid, res_1_mat),]
res_2_match = grid[whichmindist_func(grid, res_2_mat),]

#set up the kernal function -------------------
value_df = sdf(0, forage)
max_kern = max_kern_func(res_1_match, res_2_match)

#evaluate landscape -------------------
res_val = res_value_function(res_1_match, study_grid_matrix,ac)
max_val = landscape_value_sum(res_1_match, res_2_match, res_val)
cell_value = cellvalue(res_1_match, res_2_match, res_val)

#plot results -------------------
Col =  rbPal(50)[as.numeric(cut(cell_value$val,breaks = 50))]
ras_value = raster_value_func(cell_value, grid)
plot(ras_value, interpolate = T)
points(res_1_match, cex = 1.7, pch = 20, col = 'black')
points(res_2_match, cex = 1.7, pch = 20, col = 'blue')

#least cost pathing -------------------
pass = passage_func(res_1_match, res_2_match, cell_value, ras_value, perm)
plot(pass, col = cols_out, interpolate = T)

#Pareto optimisation, this can be performed on the pass, or ras_value layers -------------------
#you must click a location on the Plot section to make your decision and click finish in the top right
par_out = pareto_func(ras_value ,conf) #read above
plot(ras_value)
plot(par_out, col = 'red', legend = FALSE, add = TRUE)





























