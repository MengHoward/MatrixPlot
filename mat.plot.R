require(HiTC)
require(GenomicRanges)
require(Matrix)

chrom.name <- function(chrom_index){
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
}

chrom.name.factor <- function(chrom_index_list){
  chrom_index_list = as.vector(chrom_index_list)
  chrom_name_list = rep("chrN",length(chrom_index_list))
  for(index in c(1:length(chrom_index_list))){
    chrom_index = chrom_index_list[index]
    if(chrom_index < 23){
      chrom_name = paste0("chr",chrom_index)
    }else if(chrom_index == 23){
      chrom_name = "chrX"
    }else if(chrom_index == 24){
      chrom_name = "chrY"
    }else if(chrom_index == 25){
      chrom_name = "rDNA"
    }else{
      chrom_name = NA
    }
    chrom_name_list[index] = chrom_name
  }
  return(chrom_name_list)
}

chrom.length <- function(chrom_index,GenomeVersion="hg19"){
  # 根据chr_index返回chr_len
  HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
             135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
             51304566,155270560,59373566,42999)
  
  HG38_LEN=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
             135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
             50818468,156040895,57227415,42999)
  if(GenomeVersion=="hg19"){
    return(HG19_LEN[chrom_index])
  }else if(GenomeVersion=="hg38"){
    return(HG38_LEN[chrom_index])
  }
}

make_HTCexp_obj <- function(dt, chr_index_1,chr_index_2,resolution){
  # 构造HTCexp对象
  # 参数：
  #   dt: Hi-C连接矩阵(full)
  #   chr: 矩阵是哪条染色体上的
  #   chr.length: 染色体长度
  #   resolution: Hi-C分辨率
  # 返回：
  #   由连接矩阵生成的HTCexp对象
  
  chr_name_1 <- paste0("chr",chr_index_1)
  chr_name_2 <- paste0("chr",chr_index_2)
  
  gr.obj_x <- GRanges(
    seqnames = chr_name_2,
    ranges = IRanges(seq(0, chrom.length(chr_index_2), resolution),
                     seq(0, chrom.length(chr_index_2), resolution) + resolution),
    name = seq(0, chrom.length(chr_index_2), resolution)
  )
  
  gr.obj_y <- GRanges(
    seqnames = chr_name_1,
    ranges = IRanges(seq(0, chrom.length(chr_index_1), resolution),
                     seq(0, chrom.length(chr_index_1), resolution) + resolution),
    name = seq(0, chrom.length(chr_index_1), resolution)
  )
  
  htc.mat <- Matrix(as.matrix(dt))
  rownames(htc.mat) <- seq(0, chrom.length(chr_index_1), resolution)
  colnames(htc.mat) <- seq(0, chrom.length(chr_index_2), resolution)
  htc.obj <- new("HTCexp", htc.mat, gr.obj_x, gr.obj_y)
  return(htc.obj)
}


make.htc_exp.obj <- function(dt,chrom_index,chrom_start=0,chrom_end=0,binsize=1e6){
  # 构造同1条染色体的HTCexp对象
  # 参数：
  #   dt: Hi-C连接矩阵(full) n*n
  #   chrom_index: 矩阵是哪条染色体上的 1..22.23.24
  #   chrom_start: 矩阵的起始位置 默认 0
  #   chrom_end:  矩阵的终止位置 默认 chrom_len
  #   binsize: Hi-C matrix binsize
  # 返回：
  #   由连接矩阵生成的HTCexp对象
  chrom_len = chrom.length(chrom_index)
  
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  if(chrom_end <= chrom_start){
    chrom_end = chrom_len
  }
  
  range_start_index = chrom_start %/% binsize
  if(range_start_index == 0){range_start_index = 1}
  
  range_end_index = chrom_end %/% binsize + 1
  range_start = range_start_index * binsize
  range_end = range_end_index * binsize 
  
  gr.obj <- GRanges(
    seqnames = chrom_name,
    ranges = IRanges(seq(range_start, range_end, binsize),
                     seq(range_start, range_end, binsize) + binsize),
    name = seq(range_start, range_end, binsize)
  )
  
  htc.mat <- Matrix(as.matrix(dt)[range_start_index:range_end_index,range_start_index:range_end_index])
  rownames(htc.mat) <- seq(range_start, range_end, binsize)
  colnames(htc.mat) <- seq(range_start, range_end, binsize)
  htc.obj <- new("HTCexp", htc.mat, gr.obj, gr.obj)
  return(htc.obj)
}


matrix.fix <- function(mat,min_bound=0,max_bound=1){
  mat <- as.matrix(mat)
  # matrix value fixation
  min_bound.quantile <- quantile(mat,min_bound,type = 1)
  max_bound.quantile <- quantile(mat,max_bound,type = 1)
  mat[mat <= min_bound.quantile] <- min_bound.quantile
  mat[mat >= max_bound.quantile] <- max_bound.quantile
  return(mat)
}

matrix.part <- function(mat,chrom_index,chrom_start=0,chrom_end=NULL,binsize=1e6){
  mat = as.matrix(mat)
  chrom_len = chrom.length(chrom_index)
  if(is.null(chrom_end)){chrom_end=chrom_len}
  
  chrom_start.index = chrom_start %/% binsize + 1
  chrom_end.index = chrom_end %/% binsize + 1
  
  return(mat[chrom_start.index:chrom_end.index,chrom_start.index:chrom_end.index])
}

matrix.part.corner <- function(mat,chrom_index,segment.row=c(0,1e6),segment.col=NULL,binsize=1e6){
  mat = as.matrix(mat)
  chrom_len = chrom.length(chrom_index)
  if (is.null(segment.col)) {
    segment.col = c(0,chrom_len)
  }
  
  segment.row.start.index = segment.row[1] %/% binsize + 1
  segment.row.end.index = segment.row[2] %/% binsize + 1
  segment.col.start.index = segment.col[1] %/% binsize + 1
  segment.col.end.index = segment.col[2] %/% binsize + 1
  
  return(mat[segment.row.start.index:segment.row.end.index,
             segment.col.start.index:segment.col.end.index]
         )
}

matrix.norm.ice <- function(mat, chrom_index, binsize = 50e3, max_iter = 50){
  # 使用ICE校正Hi-C连接矩阵
  # 参数：
  #   mat: Hi-C矩阵，n*n
  #   chr: 要校正的染色体
  #   chr.size: 染色体长度
  #   resolution: Hi-C分辨率
  #   max_iter: 校正时迭代次数
  # 返回：
  #   经过校正的Hi-C矩阵
  htc.obj <- make_HTCexp_obj(mat,chr_index_1 = chrom_index,chr_index_2 = chrom_index,resolution = binsize)
  htc.obj.ice <- normICE(htc.obj, max_iter)
  hic.chr.norm <- intdata(htc.obj.ice)
  hic.chr.norm <- as.matrix(hic.chr.norm)
  return(hic.chr.norm)
}

# plot.matrix <- function(mat,min_bound=0,max_bound=1,color_type=1){
#   # mat = hic matrix
#   # min_bound min quantile to miss data
#   # max_bound max quantile to miss data
#   
#   mat <- as.matrix(mat)
#   
#   #matrix info calculate
#   row_num <- dim(mat)[1]
#   col_num <- dim(mat)[2]
#   
#   #matrix rects' coordinate
#   x1 <- rep(c(0:(col_num-1)),each=row_num)
#   x2 <- x1 + 1
#   y1 <- rep(c(0:(-row_num+1)),col_num)
#   y2 <- y1 -1
#   
#   #matrix colour vector
#   if(color_type == 1){
#     #colType == 1, white-red color pattern
#     #matrix value fixation
#     mat.fix <- matrix.fix(mat,min_bound,max_bound)
#     col_vector <- as.vector(mat.fix) / max(mat.fix)
#     red_vector <- rep(1,length(col_vector))
#     green_vector <- 1 - col_vector
#     blue_vector <- 1 - col_vector
#   }else if(color_type == 2){
#     # red blue 
#   }
#   
#   # plot the final matrix
#   plot(x=c(0,col_num),y=c(-row_num,0),type="n",frame.plot = F,xaxt="n",yaxt="n",cex.main = 2,xlab="",ylab="")
#   rect(x1,y1,x2,y2,col = rgb(red_vector,green_vector,blue_vector),border = NA)
# }


plot.matrix <- function(mat,bound.min=0,bound.max=1,mat.lim=NULL,color_type=1,col.min = "red",col.max = "red",col.boundary = NULL,n_block_color="#FFFFFF"){
  # mat = hic matrix
  # min_bound min quantile to miss data
  # max_bound max quantile to miss data
  
  mat <- as.matrix(mat)
  
  #matrix info calculate
  row_num <- dim(mat)[1]
  col_num <- dim(mat)[2]
  
  #matrix rects' coordinate
  x1 <- rep(c(0:(col_num-1)),each=row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-row_num+1)),col_num)
  y2 <- y1 -1
  
  #matrix colour vector
  if(color_type == 1){
    ## 生成矩阵需要的颜色
    ###确定 matrix的上下界
    mat.quantile = quantile(as.vector(mat),prob=c(bound.min,bound.max))
    
    if(is.null(mat.lim)){
      mat.lim = mat.quantile
    }else{
      mat.lim = c(min(c(mat.quantile,mat.lim)),max(c(mat.quantile,mat.lim)))
    }
    
    if(is.null(col.boundary)){
      col.boundary = mat.lim[1]
    }
    
    if(col.boundary>mat.lim[2]){
      print("Error! col.boundary have to smaller than matrix mat.lim[2] value!")
      return(NULL)
    }else if(col.boundary<mat.lim[1]){
      print("Error! col.boundary have to larger than matrix mat.lim[1] value!")
      print(mat.lim)
      print(col.boundary)
      return(NULL)
    }
    
    ## fix value,too large or too small
    mat[mat < mat.lim[1]] = mat.lim[1]
    mat[mat > mat.lim[2]] = mat.lim[2]
    
    ## create color vector
    mat.col_alpha = rep(0,length(as.vector(mat)))
    mat.col_alpha[mat >= col.boundary] = ceiling((mat[mat >= col.boundary] - col.boundary) / (mat.lim[2] - col.boundary) * 255)
    mat.col_alpha[mat < col.boundary] = ceiling((col.boundary - mat[mat < col.boundary]) / (col.boundary - mat.lim[1]) * 255)
    
    mat.color = rep("#FFFFFF",length(as.vector(mat)))
    mat.color[mat>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.col_alpha[mat>=col.boundary],maxColorValue = 255)
    mat.color[mat< col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.col_alpha[mat< col.boundary],maxColorValue = 255)
    
    mat.color = matrix(mat.color,nrow = row_num,ncol = col_num)
    mat.color[,colSums(mat)==0] = rgb(t(col2rgb(n_block_color)),maxColorValue = 255)
    mat.color[rowSums(mat)==0,] = rgb(t(col2rgb(n_block_color)),maxColorValue = 255)
    
  }else if(color_type == 2){
    # input a col_list 
  }
  
  # plot the final matrix
  plot(x=c(0,col_num),y=c(-row_num,0),type="n",frame.plot = F,xaxt="n",yaxt="n",cex.main = 2,xlab="",ylab="")
  rect(x1,y1,x2,y2,col = mat.color,border = NA)
}


plot.track <- function(df.peak,chrom_index,chrom_len,chrom_start=0,chrom_end=chrom_len,color_pattern = FALSE,track_color = "red",track_height = 10,axis_x_len = 50e6,xaxt = TRUE,yaxt = TRUE,xcex = 2,ycex = 2,xylwd = 5,rect_border=NA){
  
  # 参数列表
  ## df.peak dataframe 需要输入的bed narrow/broad peak 文件，格式具体参考 https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  ## chrom_index int 染色体编号，1..24，23 means X，24 means Y
  ## chrom_len int 染色体的长度
  ## chrom_start int 需要画图的区域的起始位置 default=0
  ## chrom_end int 需要画图的区域的终止位置 default=chrom_len
  ## color_pattern logical TRUE则根据bed value绘制256度灰，FALSE则绘制统一颜色; default=TRUE
  ## track_color color 当color_pattern==FALSE的时候，统一trick的颜色 default="red"
  ## track_height int 图形高度 default=10
  ## axis_x_len int x轴标记的间隔长度 default=1e6
  ## xaxt logical 是否画x轴 default=TRUE
  ## yaxt logical 是否画y轴 default=TRUE
  ## xcex x轴字体放大倍数 default=2
  ## ycex y轴字体放大倍数 default=2
  ## xylwd xy轴轴线的粗细 default=5
  
  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  # 只画某1条染色体，筛选数据
  bed_table = df.peak[df.peak[,1]==chrom_name,]
  filter_vector = (bed_table[,2] >= chrom_start & bed_table[,2] <= chrom_end) & (bed_table[,3] >= chrom_start & bed_table[,3] <= chrom_end)
  bed_table = bed_table[filter_vector,]
  
  # 只使用有用的数据信息
  track_start = bed_table[,2]
  track_end = bed_table[,3]
  track_value = bed_table[,5]
  
  # 选择color模式
  ## 如果color_pattern 为 TRUE，那么是根据track_value画256色灰度
  ## 如果color_pattern 为 FALSE，那么统一颜色，使用track_color参数
  if(color_pattern){
    ## use 255 gray pattern
    red_vector = floor(track_value / 1000 * 255)
    red_vector = 255 - red_vector
    green_vector = red_vector
    blue_vector = red_vector
    color_vector = rgb(red_vector,green_vector,blue_vector,maxColorValue = 255)
  }else{
    ## same color
    color_vector = rep(track_color,nrow(bed_table))
  }
  
  # rect 坐标
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  track_y1 <- rep((-1)*(track_height %/% 2),nrow(bed_table))
  track_y2 <- rep(track_height %/% 2,nrow(bed_table))
  
  # 绘图区域
  ## 生成画布
  plot(x=c(chrom_start,chrom_end),y=c((-1)*(track_height %/% 2),(track_height %/% 2)),
       type="n",frame.plot = F,cex.axis=2,cex.lab=2,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  ## 画矩形
  rect(track_x1,track_y1,track_x2,track_y2,col = color_vector,border = rect_border)
  
  ## 画y轴
  if(yaxt){
    axis_y_at = c(track_height %/% 2,-track_height %/% 2)
    axis(side=2,at=axis_y_at,cex.axis=ycex,lwd.ticks = xylwd,lwd = xylwd,labels = F)
  }
  
  ## 画x轴
  if(xaxt){
    if(axis_x_len >= 1e6){
      axis_x_label_size = 1e6
    }else if(axis_x_len >= 1e3){
      axis_x_label_size = 1e3
    }else{
      axis_x_label_size = 1
    }
    axis_x_at = seq(from=chrom_start,to=(chrom_end %/% axis_x_len + 1)*axis_x_len,by = axis_x_len)
    axis_x_label = axis_x_at %/% axis_x_label_size
    axis(side=1,at=axis_x_at,labels = axis_x_label,cex.axis=xcex,lwd = xylwd,lwd.ticks = xylwd)
  }
}


plot.peak <- function(df.peak,
                      chrom_index,chrom_len,chrom_start=0,chrom_end=chrom_len,
                      color_pattern = T,track_color = "red",track_pattern=F,track_height = 10,track.frame=F,
                      axis_x_len = 50e6,xaxt = TRUE,yaxt = TRUE,xcex = 2,ycex = 2,xylwd = 5,peak.lwd=1,rect_border=NA,rect_border.lwd=1){
  
  # 参数列表
  ## df.peak dataframe 需要输入的bed narrow/broad peak 文件，格式具体参考 https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  ## chrom_index int 染色体编号，1..24，23 means X，24 means Y
  ## chrom_len int 染色体的长度
  ## chrom_start int 需要画图的区域的起始位置 default=0
  ## chrom_end int 需要画图的区域的终止位置 default=chrom_len
  ## color_pattern logical TRUE则根据bed value绘制256度灰，FALSE则绘制统一颜色; default=TRUE
  ## track_color color 当color_pattern==FALSE的时候，统一trick的颜色 default="red"
  ## track_height int 图形高度 default=10
  ## axis_x_len int x轴标记的间隔长度 default=1e6
  ## xaxt logical 是否画x轴 default=TRUE
  ## yaxt logical 是否画y轴 default=TRUE
  ## xcex x轴字体放大倍数 default=2
  ## ycex y轴字体放大倍数 default=2
  ## xylwd xy轴轴线的粗细 default=5
  
  #   df.peak = active_h3k4me3_track
  #   chrom_index = 1
  #   chrom_start = 10e6
  #   chrom_end = 20e6
  #   color_pattern = T
  #   track_color = "blue"
  #   track_height = 10
  #   rect_border = F
  
  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  # 只画某1条染色体，筛选数据
  bed_table = df.peak[df.peak[,1]==chrom_name,]
  filter_vector = (bed_table[,2] >= chrom_start & bed_table[,2] <= chrom_end) & (bed_table[,3] >= chrom_start & bed_table[,3] <= chrom_end)
  bed_table = bed_table[filter_vector,]
  
  # 只使用有用的数据信息
  track_start = bed_table[,2]
  track_end = bed_table[,3]
  track_value = bed_table[,5]
  
  # 选择color模式
  ## 如果color_pattern 为 TRUE，那么是根据track_value画256色灰度
  ## 如果color_pattern 为 FALSE，那么统一颜色，使用track_color参数
  if(is.null(track_color)){
    track_color = rgb(0,0,0,maxColorValue = 255,alpha = 1)
  }
  
  if(color_pattern){
    ## use 255 gray pattern
    alpha_vector = floor(track_value / 600 * 255) # 这里设置600纯是为了好看
    alpha_vector[alpha_vector>255] = 255
    color_vector = rgb(t(col2rgb(track_color)),alpha = alpha_vector,maxColorValue = 255)
    
  }else{
    ## same color
    color_vector = rep(track_color,nrow(bed_table))
  }
  
  # rect 坐标
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  track_y1 <- rep((-1)*(track_height %/% 2),nrow(bed_table))
  track_y2 <- rep(track_height %/% 2,nrow(bed_table))
  
  # 绘图区域
  ## if track_pattern == TRUE, plot rect track
  if(track_pattern){
    plot(x=c(chrom_start,chrom_end),y=c((-1)*(track_height %/% 2),(track_height %/% 2)),
         type="n",frame.plot = track.frame,cex.axis=2,cex.lab=2,ylab = "",xlab = "",xaxt = "n",yaxt="n")
    rect(track_x1,track_y1,track_x2,track_y2,col = color_vector,border = rect_border,lwd = rect_border.lwd)
  }else{
    plot(x=track_start,y=track_value,type="h",xlim=c(chrom_start,chrom_end),ylim=c(0,1000),lwd=peak.lwd,col=track_color,frame.plot=track.frame,xaxt = "n",yaxt="n",ylab = "",xlab = "")
  }
  
  ## 画y轴
  if(yaxt){
    if(track_pattern){
      axis_y_at = c(track_height %/% 2,0,-track_height %/% 2)  
    }else{
      axis_y_at = c(0,500,1000)  
    }
    
    axis(side=2,at=axis_y_at,cex.axis=ycex,lwd.ticks = xylwd,lwd = xylwd,labels = F)
  }
  
  ## 画x轴
  if(xaxt){
    if(axis_x_len >= 1e6){
      axis_x_label_size = 1e6
    }else if(axis_x_len >= 1e3){
      axis_x_label_size = 1e3
    }else{
      axis_x_label_size = 1
    }
    axis_x_at = seq(from=chrom_start,to=(chrom_end %/% axis_x_len + 1)*axis_x_len,by = axis_x_len)
    axis_x_label = axis_x_at %/% axis_x_label_size
    axis(side=1,at=axis_x_at,labels = axis_x_label,cex.axis=xcex,lwd = xylwd,lwd.ticks = xylwd)
  }
}

gene.density <- function(gene.df,chrom_index,chrom_start,chrom_end){
  # 计算区域内的gene density
  # gene.df: 给出的gene bed格式
  # chrom_index: 染色体编号1..22,23,24
  # chrom_start: 区域的起始位置
  # chrom_end: 区域的终止位置
  
  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  filter_vector <- (gene.df$chrom_name == chrom_name) & (gene.df$start >= chrom_start) & (gene.df$start <= chrom_end)
  return(length(which(filter_vector==T)) / (chrom_end - chrom_start))
}


overlap <- function(table_1,table_2,level_vetor=c("chr1")){
  # overlap 是针对排序过的table_1,table_2进行取overlap
  # 返回是1个overlap的data.frame
  overlap_table = data.frame(row.names = c("chrom_name","start","end","name","start.1","end.1","name.1","start.2","end.2","name.2"))
  
  for(table.level in level_vetor){
    # 选择子集
    bed_table_1 = table_1[table_1[,1]==table.level,]
    bed_table_2 = table_2[table_2[,1]==table.level,]
    
    # 开始循环
    index_1 = 1
    index_2 = 1
    overlap_index = 1
    
    while(index_1 <= nrow(bed_table_1) & index_2 <= nrow(bed_table_2)){
      region_1.start = bed_table_1[index_1,2]
      region_1.end = bed_table_1[index_1,3]
      region_1.name = bed_table_1[index_1,4]
      
      region_2.start = bed_table_2[index_2,2]
      region_2.end = bed_table_2[index_2,3]
      region_2.name = bed_table_2[index_2,4]
      
      if(region_1.end < region_2.start){
        # if region 1 is less than the region 2
        index_1 = index_1 + 1
      }else if(region_1.start > region_2.end){
        # if region 2 is less than the region 1
        index_2 = index_2 + 1
      }else{
        # must have overlap region
        overlap.set = sort(c(region_1.start,region_1.end,region_2.start,region_2.end))
        overlap.start = overlap.set[2]
        overlap.end = overlap.set[3]
        overlap_row = data.frame(chrom_name = table.level,
                                 start = overlap.start,
                                 end = overlap.end,
                                 name = sprintf("overlap-%d",overlap_index),
                                 start.1 = region_1.start,
                                 end.1 = region_1.end,
                                 name.1 = region_1.name,
                                 start.2 = region_2.start,
                                 end.2 = region_2.end,
                                 name.2 = region_2.name
        )
        overlap_table = rbind(overlap_table,overlap_row)
        overlap_index = overlap_index + 1
        
        if(region_1.end < region_2.end){
          index_1 = index_1 + 1
        }else{
          index_2 = index_2 + 1
        }
      }
    }
  }
  return(overlap_table)
}


hg19_g_band = read.table("~/menghw_HD/reference/hg19_G-band.txt",header = F,sep = "\t")
plot.chromosome <- function(chrom_index,chrom_start=0,chrom_end=250e6,track_height=8,border.lwd = 1,genome_version="hg19"){
  # initialize
  # chrom_index = 1
  # chrom_start = 0
  # chrom_end = chrom.length(chrom_index)
  df.peak = hg19_g_band

  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  # 只画某1条染色体，筛选数据
  bed_table = df.peak[df.peak[,1]==chrom_name,]
  filter_vector = (bed_table[,2] >= chrom_start & bed_table[,2] <= chrom_end) & (bed_table[,3] >= chrom_start & bed_table[,3] <= chrom_end)
  bed_table = bed_table[filter_vector,]
  
  # 只使用有用的数据信息
  track_start = bed_table[,2]
  track_end = bed_table[,3]
  track_value = bed_table[,5]
  track_info = bed_table[,6]

  ## use 255 gray pattern
  alpha_vector = floor(track_value / 600 * 255) # 这里设置600纯是为了好看
  alpha_vector[alpha_vector>255] = 255
  color_vector = rgb(t(col2rgb("black")),alpha = alpha_vector,maxColorValue = 255)
    
  # rect 坐标
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  track_y1 <- rep((-1)*(track_height %/% 2),nrow(bed_table))
  track_y2 <- rep(track_height %/% 2,nrow(bed_table))
  
  track_y1[track_info=="acen"] <- (-1)*(track_height %/% 2)/2
  track_y2[track_info=="acen"] <- (track_height %/% 2)/2
  
  plot(x=c(chrom_start,chrom_end),y=c((-1)*(track_height %/% 2),(track_height %/% 2)),type="n",frame.plot = F,cex.axis=2,cex.lab=2,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  rect(track_x1,track_y1,track_x2,track_y2,col = color_vector,border = T,lwd = border.lwd)
}



complement <- function(table_1,table_2,level_vetor=c("chr1")){
  # complement 是针对排序的table，从table1中取得table2的补集
  # 返回是1个data.frame
  complement_table = data.frame(row.names = c("chrom_name","start","end","name","value","region_length"))
  
  for(table.level in level_vetor){
    # 选择子集
    bed_table_1 = table_1[table_1[,1]==table.level,]
    bed_table_2 = table_2[table_2[,1]==table.level,]
    
    # 开始循环
    index_1 = 1
    index_2 = 1
    complement_index = 1
    
    region_1.read_state = T
    region_2.read_state = T
    
    while(index_1 <= nrow(bed_table_1) & index_2 <= nrow(bed_table_2)){
      
      if(region_1.read_state == T){
        region_1.start = bed_table_1[index_1,2]
        region_1.end = bed_table_1[index_1,3]
        region_1.read_state = F
      }
      
      if(region_2.read_state == T){
        region_2.start = bed_table_2[index_2,2]
        region_2.end = bed_table_2[index_2,3]  
        region_2.read_state = F
      }
      
      if(region_1.end < region_2.start){
        # if region 1 is less than the region 2
        comlement_row = data.frame(chrom_name = table.level,
                                   start = region_1.start,
                                   end = region_1.end,
                                   name = sprintf("index-%d",complement_index),
                                   value = 1,
                                   region_length = region_1.end - region_1.start + 1)
        complement_table = rbind(complement_table,comlement_row)
        complement_index = complement_index + 1
        index_1 = index_1 + 1
        region_1.read_state = T
        
      }else if(region_1.start > region_2.end){
        # if region 2 is less than the region 1
        index_2 = index_2 + 1
        region_2.read_state = T
        
      }else{
        # must have overlap region
        if(region_1.end >= region_2.end){
          # region 2 have been read out
          if(region_1.start >= region_2.start){
            region_1.start = region_2.end + 1
          }else if(region_1.start < region_2.start){
            comlement_row = data.frame(chrom_name = table.level,
                                     start = region_1.start,
                                     end = region_2.start-1,
                                     name = sprintf("index-%d",complement_index),
                                     value = 1,
                                     region_length = region_2.start - region_1.start)
            complement_table = rbind(complement_table,comlement_row)
            complement_index = complement_index + 1
            region_1.start = region_2.end + 1
          }
          index_2 = index_2 + 1
          region_2.read_state = T
          
        }else if(region_1.end < region_2.end){
          # region 1 have been read out
          if(region_1.start <= region_2.start){
            comlement_row = data.frame(chrom_name = table.level,
                                       start = region_1.start,
                                       end = region_2.start-1,
                                       name = sprintf("index-%d",complement_index),
                                       value = 1,
                                       region_length = region_2.start - region_1.start)
            complement_table = rbind(complement_table,comlement_row)
            complement_index = complement_index + 1
            region_2.start = region_1.end + 1
          }else if(region_1.start > region_2.start){
            region_2.start = region_1.end + 1
          }
          index_1 = index_1 + 1
          region_1.read_state = T
        }
      }
    }
    
  }
  return(complement_table)
}

# hg19_N_pattern = read.table("~/menghw_HD/reference/hg19_N_pattern.bed",header = F,sep = "\t")
# hg19_length = read.table("~/menghw_HD/reference/hg19_length.bed",header = F,sep = "\t")
# hg19_not_n_table = complement(table_1 = hg19_length,table_2 = hg19_N_pattern,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# write.table(hg19_not_n_table,file = "~/menghw_HD/reference/hg19_not_N_pattern.bed",quote = F,col.names = T,row.names = F,sep = "\t")

# save.image(file = "~/menghw_HD/R_code/my_function/MyPlot_v03.RData")



##########################################################################################
# 第4次修改，增加matrix to WashU pairwise interaction file
# save.image(file = "~/menghw_HD/R_code/my_function/MyPlot_v04.RData")
# 
# 
# 
# chrom.length <- function(chrom_index,GenomeVersion="hg19"){
#   # 根据chr_index返回chr_len
#   HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
#              135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
#              51304566,155270560,59373566,42999)
# 
#   HG38_LEN=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
#              135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
#              50818468,156040895,57227415,42999)
#   if(GenomeVersion=="hg19"){
#     return(HG19_LEN[chrom_index])
#   }else if(GenomeVersion=="hg38"){
#     return(HG38_LEN[chrom_index])
#   }
# }
# 
# 
# 
# convert_matrix2inter.cis <- function(hic_matrix,chrom_index,chrom_start,chrom_end,binsize){
# 
#   if(chrom_index <=22){
#     chrom_name = paste0("chr",chrom_index)
#   }else if(chrom_index == 23){
#     chrom_name = "chrX"
#   }else if(chrom_index == 24){
#     chrom_name = "chrY"
#   }else if(chrom_index == 25){
#     chrom_name = "rDNA"
#   }
#   
#   chrom_loci_index.start = chrom_start %/% binsize + 1
#   chrom_loci_index.end = chrom_end %/% binsize + 1
# 
#   # hic_matrix.table = data.frame(row.names = c("chrom_name.1","start.1","end.1","chrom_name.2","start.2","end.2","value","value_log2"))
#   
#   vector.len = (chrom_loci_index.end - chrom_loci_index.start + 1 + 1) * (chrom_loci_index.end - chrom_loci_index.start + 1) / 2
#   
#   # chrom_name.1.vector = NULL
#   chrom_name.1.vector = rep(chrom_name,vector.len)
#   start.1.vector= NULL
#   end.1.vector= NULL
#   # chrom_name.2.vector= NULL
#   chrom_name.2.vector= rep(chrom_name,vector.len)
#   start.2.vector= NULL
#   end.2.vector= NULL
#   value.vector= NULL
#   value_log2.vector= NULL
#   
#   for(i in c(chrom_loci_index.start:chrom_loci_index.end)){
#     loci.start.1 = rep((i-1) * binsize,chrom_loci_index.end - i + 1)
#     loci.start.2 = (c(i:chrom_loci_index.end) - 1) * binsize
#     
#     loci.end.1 = loci.start.1 + binsize - 1
#     loci.end.2 = loci.start.2 + binsize - 1
#     
#     start.1.vector = c(start.1.vector,loci.start.1)
#     end.1.vector = c(end.1.vector,loci.end.1)
#   
#     start.2.vector = c(start.2.vector,loci.start.2)
#     end.2.vector = c(end.2.vector,loci.end.2)
#   }
#   
#   hic_part = as.vector(hic_matrix[lower.tri(as.matrix(hic_matrix[chrom_loci_index.start:chrom_loci_index.end,chrom_loci_index.start:chrom_loci_index.end]),diag = T)])
#   
#   hic_matrix.part = hic_matrix[chrom_loci_index.start:chrom_loci_index.end,chrom_loci_index.start:chrom_loci_index.end]
#   value.vector = hic_matrix.part[lower.tri(as.matrix(hic_matrix.part),diag = T)]
#   value_log2.vector = log2(value.vector + 1)
#   
#   hic_matrix.table = data.frame(chrom_name.1 = chrom_name.1.vector,
#                          start.1 = start.1.vector,
#                          end.1 = end.1.vector,
#                          chrom_name.2 = chrom_name.2.vector,
#                          start.2 = start.2.vector,
#                          end.2 = end.2.vector,
#                          value = value.vector,
#                          value_log2 = value_log2.vector)
# 
#   return(hic_matrix.table)
# }
# 
# hic_matrix = read.table("~/menghw_HD/our_data/Hela-genome-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_20000/chr_20_20000_MAPQ20.txt",header = F,sep = ",")
# hic_inter_table = convert_matrix2inter.cis(hic_matrix,20,10e6,20e6,20e3)
# 
# 
# # 
# 
# convert_inter2file.cis <- function(hic_inter_table,log_scale=F){
#   log_scale = F
#   WashU_pairwise.fmt = "%s:%d,%d\t%s:%d,%d\t%.4f"
#   write_vector = NULL
#   
#   for(i in c(1:nrow(hic_inter_table))){
#     print(i)
#     inter_table.row = hic_inter_table[i,]
#     inter_table.str = NULL
#     if(log_scale){
#       inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,8])
#     }else{
#       inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,7])
#     }
#     write_vector = c(write_vector,inter_table.str)
#   }
#   
#   write(x = write_vector,file = "~/menghw_HD/test.washu")
# }
# 
# 
# 
# covert_matrix2file.cis <- function(hic_matrix,chrom_index,chrom_start,chrom_end,binsize,log_scale=F,file_path){
#   
#   if(chrom_index <=22){
#     chrom_name = paste0("chr",chrom_index)
#   }else if(chrom_index == 23){
#     chrom_name = "chrX"
#   }else if(chrom_index == 24){
#     chrom_name = "chrY"
#   }else if(chrom_index == 25){
#     chrom_name = "rDNA"
#   }
#   
#   chrom_loci_index.start = chrom_start %/% binsize
#   chrom_loci_index.end = chrom_end %/% binsize
#   
#   WashU_pairwise.fmt = "%s:%d,%d\t%s:%d,%d\t%.4f"
#   write_vector = NULL
#   hic_matrix.table = data.frame(row.names = c("chrom_name.1","start.1","end.1","chrom_name.2","start.2","end.2","value","value_log2"))
#   
#   for(i in c(chrom_loci_index.start+1:chrom_loci_index.end+1)){
#     for(j in c(i:chrom_loci_index.end+1)){
#       loci.start.1 = i * binsize
#       loci.start.2 = j * binsize
#       loci.end.1 = loci.start.1 + binsize - 1
#       loci.end.2 = loci.start.2 + binsize - 1
#       
#       loci.interaction.raw = hic_matrix[i,j]
#       
#       if(loci.interaction.raw == 0){
#         loci.interaction.log2 = 0
#       }else{
#         loci.interaction.log2 = log2(loci.interaction.raw)
#       }
#       
#       
#       write_vector = c(write_vector,sprintf(WashU_pairwise.fmt,chrom_name,loci.start.1,loci.end.1,chrom_name,loci.start.2,loci.end.2,loci.interaction))
#     }
#   }
#   
#   write(x=write_vector,file = file_path)
# }

plot.barplot <- function(bed_table_all,chrom_index,chrom_start=0,chrom_end=0,ylim=NULL,track_color = "red"){
  
  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  # 只画某1条染色体，筛选数据
  bed_table = bed_table_all[bed_table_all[,1]==chrom_name,]
  filter_vector = (bed_table[,2] >= chrom_start & bed_table[,2] <= chrom_end) & (bed_table[,3] >= chrom_start & bed_table[,3] <= chrom_end)
  bed_table = bed_table[filter_vector,]
  
  # 只使用有用的数据信息
  track_start = bed_table[,2]
  track_end = bed_table[,3]
  track_value = bed_table[,4]
  
  # 选择color
  color_vector = rep(track_color,nrow(bed_table))
  
  # rect 坐标
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  
  if(is.null(ylim)){
    ylim = c(0,max(track_value))
  }
  
  track_y1 <- rep(ylim[1],nrow(bed_table))
  track_y2 <- track_value
  
  # 绘图区域
  ## if track_pattern == TRUE, plot rect track
  plot(x=c(chrom_start,chrom_end),y=ylim,type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  rect(track_x1,track_y1,track_x2,track_y2,col = color_vector,border = F,lwd = 0)
  
}


# plot.TAD <- function(input.matrix,color="red",ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T){
#   
#   # 参数说明
#   ## input.matrix: 输入矩阵
#   ## color="red": 热图颜色，支持字符或者是RGB参数
#   ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
#   ## minBound=0: 使用分位数来修正最小值
#   ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
#   ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分
#   
#   mat.nrow = dim(input.matrix)[1]
#   # 生成矩阵需要的颜色
#   mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = T)])
#   mat.quantile = quantile(mat.lower.tri,prob=c(minBound,maxBound))
#   ## fix value,too large or too small
#   mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
#   mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
#   ## create color vector 
#   mat.color.alpha = ceiling(mat.lower.tri / mat.quantile[2] * 255)
#   mat.color = rgb(t(col2rgb(color)),alpha = mat.color.alpha,maxColorValue = 255)
#   
#   # 根据矩阵的行列，生成x1,y1的坐标
#   x1.part = c(0:(mat.nrow -1))
#   y1.part = c(0:(mat.nrow -1))
#   x1 = x1.part
#   y1 = y1.part
#   for(i in c(1:(mat.nrow-1))){
#     x1.part = x1.part[-length(x1.part)] + 2
#     y1.part = y1.part[-length(y1.part)]
#     
#     x1 = c(x1,x1.part)
#     y1 = c(y1,y1.part)
#   }
#   
#   # 根据x1 y1 生成剩余坐标
#   x2 = x1 + 1
#   x3 = x2 + 1
#   x4 = x2
#   y2 = y1 + 1
#   y3 = y1
#   y4 = y1 - 1
#   # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
#   NA_vector = rep(NA,length(x1))
#   # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
#   x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
#   y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
#   x = as.vector(t(x_matrix))
#   y = as.vector(t(y_matrix))
#   
#   # 绘图部分
#   ## 如果mat.upper 为F 则画倒置的图像
#   if(! mat.upper){ y = -y }
#   ## 设置默认ylim
#   if(is.null(ylim)){ ylim =c(min(y,na.rm = T),max(y,na.rm = T)) }
#   ## 生成画布
#   plot(c(min(x,na.rm = T),max(x,na.rm = T)),c(min(y,na.rm = T),max(y,na.rm = T)),ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
#   ## 画三角矩阵
#   polygon(x,y,col = mat.color,border = F)
# }



overlap.gene <- function(region_table,gene_table,chrom_name_list=c(paste("chr",c(1:22),sep = ""),"chrX")){
  
  overlap.gene_table = NULL
  for(chrom_name in chrom_name_list){
    print(chrom_name)
    # 选择子集
    chrom_overlap.gene_table = NULL
    chrom_region_table = region_table[region_table[,1]==chrom_name,]
    chrom_gene_table = gene_table[gene_table[,1]==chrom_name,]
    for(chrom_region_index in c(1:nrow(chrom_region_table))){
      region_start = chrom_region_table[chrom_region_index,2]
      region_end = chrom_region_table[chrom_region_index,3]
      region_filter.start = (chrom_gene_table[,2] >= region_start) & (chrom_gene_table[,2] <= region_end)
      # region_filter.end = (chrom_gene_table[,3] >= region_start) & (chrom_gene_table[,3] <= region_end)
      # region_filter = region_filter.start | region_filter.end
      region_filter = region_filter.start
      region_gene_table = chrom_gene_table[region_filter,]
      chrom_overlap.gene_table = rbind(chrom_overlap.gene_table,region_gene_table)
    }
    overlap.gene_table = rbind(overlap.gene_table,chrom_overlap.gene_table)  
  }
  
  return(overlap.gene_table)
}


# plot.TAD <- function(input.matrix,color="red",ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T,mat.part=T){
#   
#   # 参数说明
#   ## input.matrix: 输入矩阵
#   ## color="red": 热图颜色，支持字符或者是RGB参数
#   ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
#   ## minBound=0: 使用分位数来修正最小值
#   ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
#   ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分; 只有在mat.part=T时有效
#   ## mat.part=T: 如果为TRUE则画三角，如果为FLASE则画倒置的matrix
# 
#   input.matrix = as.matrix(input.matrix)
#   mat.nrow = dim(input.matrix)[1]
#   
#   ## 生成矩阵需要的颜色
#   mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = F)])
#   mat.quantile = quantile(as.vector(input.matrix[lower.tri(input.matrix,diag = T)]),prob=c(minBound,maxBound))
#   ## fix value,too large or too small
#   mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
#   mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
#   mat.diagnal.value = diag(input.matrix)
#   mat.diagnal.value[mat.diagnal.value < mat.quantile[1]] = mat.quantile[1]
#   mat.diagnal.value[mat.diagnal.value > mat.quantile[2]] = mat.quantile[2]
#   
#   ## create color vector
#   mat.lower.tri.col_alpha = ceiling(mat.lower.tri / mat.quantile[2] * 255)
#   mat.diagnal.col_alpha = ceiling(mat.diagnal.value / mat.quantile[2] * 255)
#   
#   # 根据矩阵的行列，生成x1,y1的坐标
#   x1.part = c(1:(mat.nrow -1))
#   y1.part = c(1:(mat.nrow -1))
#   x1 = x1.part
#   y1 = y1.part
#   for(i in c(1:(mat.nrow-1))){
#     x1.part = x1.part[-length(x1.part)] + 2
#     y1.part = y1.part[-length(y1.part)]
#     
#     x1 = c(x1,x1.part)
#     y1 = c(y1,y1.part)
#   }
#   
#   # 对角线的x1
#   x1.diagonal = c(0:(mat.nrow-1)) * 2
#   y1.diagonal = rep(0,mat.nrow)
#   
#   if(mat.part==T){
#     # plot matrix upper OR lower only.
#     x1 = c(x1,x1.diagonal)
#     y1 = c(y1,y1.diagonal)  
#     # 颜色向量也不同
#     mat.color = rgb(t(col2rgb(color)),alpha = c(mat.lower.tri.col_alpha,mat.diagnal.col_alpha),maxColorValue = 255)
#   }else if(mat.part==F){
#     x1 = c(x1,x1.diagonal,x1)
#     y1 = c(y1,y1.diagonal,-y1)  
#     mat.color = rgb(t(col2rgb(color)),alpha = c(mat.lower.tri.col_alpha,mat.diagnal.col_alpha,mat.lower.tri.col_alpha),maxColorValue = 255)
#   }
#   
#   
#   # 根据x1 y1 生成剩余坐标
#   x2 = x1 + 1
#   x3 = x2 + 1
#   x4 = x2
#   y2 = y1 + 1
#   y3 = y1
#   y4 = y1 - 1
#   # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
#   NA_vector = rep(NA,length(x1))
#   # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
#   x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
#   y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
#   x = as.vector(t(x_matrix))
#   y = as.vector(t(y_matrix))
#   
#   # 绘图部分
#   ## 如果mat.upper 为F 则画倒置的图像
#   ## 设置画布大小
#   x.point = c(0,mat.nrow*2)
#   if(mat.part){
#     if(! mat.upper){ y = -y ; y.point = c(-mat.nrow,1)}else{ y.point = c(-1,mat.nrow) }
#   }else{ 
#     y.point = c(-mat.nrow,mat.nrow)
#   }
#   ## 设置ylim
#   if(is.null(ylim)){ ylim = y.point }
#   
#   ## 生成画布
#   plot(x.point,y.point,ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
#   ## 画三角矩阵
#   polygon(x,y,col = mat.color,border = F)
# }


# save.image(file="~/menghw_HD/R_code/my_function/MyPlot_v04.RData")

plot.TAD <- function(input.matrix,col.min = "red",col.max = "red",col.boundary = 0,ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T,mat.part=T){
  
  # 参数说明
  ## input.matrix: 输入矩阵
  ## col.min="red": 热图颜色，支持字符或者是RGB参数
  ## col.max="red": 热图颜色，支持字符或者是RGB参数
  ## col.boundary=0: 两种颜色的分界线，也即白色对应的值
  ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
  ## minBound=0: 使用分位数来修正最小值
  ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
  ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分; 只有在mat.part=T时有效
  ## mat.part=T: 如果为TRUE则画三角，如果为FLASE则画倒置的matrix
  input.matrix = as.matrix(input.matrix)
  mat.nrow = dim(input.matrix)[1]
  
  ## color type 1：single color
  ## color type 2：blue white red black
  
  ## 生成矩阵需要的颜色
  mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = F)])
  mat.quantile = quantile(as.vector(input.matrix[lower.tri(input.matrix,diag = T)]),prob=c(minBound,maxBound))
  
  if(col.boundary>mat.quantile[2] | col.boundary<mat.quantile[1]){
    print("Error! col.boundary have to smaller than matrix maxBound value!")
    return(NULL)
  }
  
  ## fix value,too large or too small
  mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
  mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
  mat.diagnal.value = diag(input.matrix)
  mat.diagnal.value[mat.diagnal.value < mat.quantile[1]] = mat.quantile[1]
  mat.diagnal.value[mat.diagnal.value > mat.quantile[2]] = mat.quantile[2]
  
  ## create color vector
  mat.lower.tri.col_alpha = rep(0,length(mat.lower.tri))
  mat.diagnal.col_alpha = rep(0,length(mat.diagnal.value))
  
  mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary] = ceiling((mat.lower.tri[mat.lower.tri>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
  mat.lower.tri.col_alpha[mat.lower.tri<col.boundary] = ceiling((col.boundary - mat.lower.tri[mat.lower.tri<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
  mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary] = ceiling((mat.diagnal.value[mat.diagnal.value>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
  mat.diagnal.col_alpha[mat.diagnal.value<col.boundary] = ceiling((col.boundary - mat.diagnal.value[mat.diagnal.value<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
  
  mat.lower.tri.color = rep("##FFFFFF",length(mat.lower.tri))
  mat.diagnal.color = rep("##FFFFFF",length(mat.diagnal.value))
  
  mat.lower.tri.color[mat.lower.tri>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary],maxColorValue = 255)
  mat.lower.tri.color[mat.lower.tri<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.lower.tri.col_alpha[mat.lower.tri<col.boundary],maxColorValue = 255)
  mat.diagnal.color[mat.diagnal.value>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary],maxColorValue = 255)
  mat.diagnal.color[mat.diagnal.value<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.diagnal.col_alpha[mat.diagnal.value<col.boundary],maxColorValue = 255)

  # 根据矩阵的行列，生成x1,y1的坐标
  x1.part = c(1:(mat.nrow -1))
  y1.part = c(1:(mat.nrow -1))
  x1 = x1.part
  y1 = y1.part
  for(i in c(1:(mat.nrow-1))){
    x1.part = x1.part[-length(x1.part)] + 2
    y1.part = y1.part[-length(y1.part)]
    
    x1 = c(x1,x1.part)
    y1 = c(y1,y1.part)
  }
  
  # 对角线的x1
  x1.diagonal = c(0:(mat.nrow-1)) * 2
  y1.diagonal = rep(0,mat.nrow)
  
  if(mat.part==T){
    # plot matrix upper OR lower only.
    x1 = c(x1,x1.diagonal)
    y1 = c(y1,y1.diagonal)  
    # 颜色向量也不同
    mat.color = c(mat.lower.tri.color,mat.diagnal.color)
  }else if(mat.part==F){
    x1 = c(x1,x1.diagonal,x1)
    y1 = c(y1,y1.diagonal,-y1)  
    mat.color = c(mat.lower.tri.color,mat.diagnal.color,mat.lower.tri.color)
  }
  
  
  # 根据x1 y1 生成剩余坐标
  x2 = x1 + 1
  x3 = x2 + 1
  x4 = x2
  y2 = y1 + 1
  y3 = y1
  y4 = y1 - 1
  # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
  NA_vector = rep(NA,length(x1))
  # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
  x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
  y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
  x = as.vector(t(x_matrix))
  y = as.vector(t(y_matrix))
  
  # 绘图部分
  ## 如果mat.upper 为F 则画倒置的图像
  ## 设置画布大小
  x.point = c(0,mat.nrow*2)
  if(mat.part){
    if(! mat.upper){ y = -y ; y.point = c(-mat.nrow,1)}else{ y.point = c(-1,mat.nrow) }
  }else{ 
    y.point = c(-mat.nrow,mat.nrow)
  }
  ## 设置ylim
  if(is.null(ylim)){ ylim = y.point }
  
  ## 生成画布
  plot(x.point,y.point,ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
  ## 画三角矩阵
  polygon(x,y,col = mat.color,border = F)
}


plot.matrix <- function(mat,bound.min=0,bound.max=1,mat.lim=NULL,color_type=1,col.min = "red",col.max = "red",col.boundary = NULL){
  # mat = hic matrix
  # min_bound min quantile to miss data
  # max_bound max quantile to miss data
  
  mat <- as.matrix(mat)
  
  #matrix info calculate
  row_num <- dim(mat)[1]
  col_num <- dim(mat)[2]
  
  #matrix rects' coordinate
  x1 <- rep(c(0:(col_num-1)),each=row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-row_num+1)),col_num)
  y2 <- y1 -1
  
  #matrix colour vector
  if(color_type == 1){
    ## 生成矩阵需要的颜色
    ###确定 matrix的上下界
    mat.quantile = quantile(as.vector(mat),prob=c(bound.min,bound.max))
    
    if(is.null(mat.lim)){
      mat.lim = mat.quantile
    }else{
      mat.lim = c(min(c(mat.quantile,mat.lim)),max(c(mat.quantile,mat.lim)))
    }
    
    if(is.null(col.boundary)){
      col.boundary = mat.lim[1]
    }
    
    if(col.boundary>mat.lim[2]){
      print("Error! col.boundary have to smaller than matrix mat.lim[2] value!")
      return(NULL)
    }else if(col.boundary<mat.lim[1]){
      print("Error! col.boundary have to larger than matrix mat.lim[1] value!")
      print(mat.lim)
      print(col.boundary)
      return(NULL)
    }
    
    ## fix value,too large or too small
    mat[mat < mat.lim[1]] = mat.lim[1]
    mat[mat > mat.lim[2]] = mat.lim[2]
    
    ## create color vector
    mat.col_alpha = rep(0,length(as.vector(mat)))
    mat.col_alpha[mat >= col.boundary] = ceiling((mat[mat >= col.boundary] - col.boundary) / (mat.lim[2] - col.boundary) * 255)
    mat.col_alpha[mat < col.boundary] = ceiling((col.boundary - mat[mat < col.boundary]) / (col.boundary - mat.lim[1]) * 255)
    
    mat.color = rep("#FFFFFF",length(as.vector(mat)))
    mat.color[mat>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.col_alpha[mat>=col.boundary],maxColorValue = 255)
    mat.color[mat< col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.col_alpha[mat< col.boundary],maxColorValue = 255)
    
  }else if(color_type == 2){
    # input a col_list 
  }
  
  # plot the final matrix
  plot(x=c(0,col_num),y=c(-row_num,0),type="n",frame.plot = F,xaxt="n",yaxt="n",cex.main = 2,xlab="",ylab="")
  rect(x1,y1,x2,y2,col = mat.color,border = NA)
}

# save.image(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

plot.segment <- function(segment.row=c(0,1),segment.col=c(0,1),negtive=T,lty=2,lwd=3,color="black"){
  x0 = c(segment.col[1],segment.col[1],segment.col[2],segment.col[2])
  x1 = c(segment.col[1],segment.col[2],segment.col[2],segment.col[1])
  
  y0 = c(segment.row[1],segment.row[2],segment.row[2],segment.row[1])
  y1 = c(segment.row[2],segment.row[2],segment.row[1],segment.row[1])
  
  if(negtive==T){
    y0 = -1 * y0
    y1 = -1 * y1
  }
  
  segments(x0,y0,x1,y1,lty=lty,lwd=lwd,col = color)
}



# save.image(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
