# example:
# hclust_out      <- hclust(as.dist(1-correlation_matrix))
# plot_clust_join_density(hclust_out)
plot_clust_join_density <- function(hclust_out) {
    
    # Let's do this a bit rigourously
    # create data frame of joining data
    df_hclust <- data.frame(x1=hclust_out$merge[,1],x2=hclust_out$merge[,2],h=hclust_out$height)
        # note that for each time two clusters are joined, this is stored
        # as a row in merge with points i,j where i is the index of one point and j the
        # the index of the other. Singletons are given negative numbers, 
        # and conglomerate clusters are given positive numbers.
        # height gives the corresponding distance of points being joints.

    # now look at increase of height    
    df_height  <- data.frame(#nr_points=seq(length(heights),1,-1),
                    nr_points=seq(1,length(df_hclust$h)),
                    heights=df_hclust$h)
    
    ################################################################################
    # 1st method to define a cutoff distance

    # In an ideal situation, there would be clearly distinct clusters in the data
    # In this case, all points first are joined to clusters, and then only at the 
    # end are the clusters joined together. At that point, the clusters are separated
    # much farther than the typical point-to-point distance, or point-to-cluster distance.
    # 
    # Thus, we can look at the distance between the last points that were merged into
    # clusters to define a cutoff distance.
    
    # get distances where single points are joined at 
    singleton_heights <- hclust_out$height[(hclust_out$merge[,1]<0)|(hclust_out$merge[,2]<0)]
    # now take max of that as a distance where "clusters" have formed
    # (or, actually, take the 98% percentile to avoid outliers)
    cutoff_distance1 <-
        sort(singleton_heights,decreasing = F)[round(length(singleton_heights)*.98)]
    cutoff_distance2 <- max(singleton_heights)
    
    # let's look at how connectedness increases
    p1<-ggplot(data=df_hclust)+
        geom_point(aes(x=x1,y=h))+
        geom_point(aes(x=x2,y=h),color='red')+
        geom_hline(yintercept = cutoff_distance1)+
        geom_hline(yintercept = cutoff_distance2)+
        ggtitle(paste0('Positive x: conglomerate joining\nNegative x: singleton joining\nHeights of lines: ',
                        round(cutoff_distance1,4),' and ', round(cutoff_distance2,4) ,'.'))+
        give_better_textsize_plot(15)+xlab('Element index')+ylab('Joining distance')
    
    ################################################################################
    # 2nd method to detect a good cutoff distance
    
    # This method tries to find the 'tipping point', where the curve goes from
    # slowly rising (short distances being added) to steeply rising (large distances
    # being added). The idea would be similar to the first method, which is that
    # the typical distance between points in clusters is smaller than the typical
    # distance between clusters of points. 
    
    # give shorthand names to heights and points
    x=df_height$nr_points
    y=df_height$heights
    
    # Extrapolate data between points
    approx_out <- approx(x=x, y=y, method='linear', n = 100)
    
    # Put the extrapolated line into convenient parameters.
    x_ = approx_out$x
    y_ = approx_out$y
    
    # calculate the derivative line of the extrapolated line
    d_y = (y_[2:length(y_)]-y_[1:(length(y_)-1)]) # just the point-to-point difference (dy)
    x_mid = (x_[2:length(y_)]+x_[1:(length(y_)-1)])/2 # the x-locations
    dy_div_dx = (y_[2:length(y_)]-y_[1:(length(y_)-1)])/(x_[2]-x_[1]) # actual derivative; dy divided by dx
    
    # calculate the 2nd derivative
    dd_y = (d_y[2:length(d_y)]-d_y[1:(length(d_y)-1)])
    x_mm = (x_mid[2:length(x_mid)]+x_mid[1:(length(x_mid)-1)])/2
    
    # determine tipping point, which should be reflected by a bump in the 2nd derivative
    # (in an ideal situation)
    # we define a bump by the cumulative sum reaching a certain percentage (e.g. 1/100) of the total area
    x_cutoff1<-x_mid[which(cumsum(dd_y)>sum(dd_y)/100)[1]]
    y_cutoff1<-y_[1+which(cumsum(dd_y)>sum(dd_y)/100)[1]]
    x_cutoff2<-x_mid[which(cumsum(dd_y)>sum(dd_y)/10)[1]]
    y_cutoff2<-y_[1+which(cumsum(dd_y)>sum(dd_y)/10)[1]]
    
    # now plot all of this in a figure
    p2<-ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))+
            geom_point(data=data.frame(x=approx_out$x, y=approx_out$y), aes(x=x,y=y,color='extrapolated'))+
            geom_point(aes(color='data'))+
            geom_vline(xintercept = x_cutoff1)+geom_hline(yintercept = y_cutoff1)+
            geom_vline(xintercept = x_cutoff2)+geom_hline(yintercept = y_cutoff2)+
            geom_text(aes(x=x_cutoff1,y=max(y),label=round(x_cutoff1,3)),hjust=1,color='red')+
            geom_text(aes(x=x_cutoff2,y=max(y)*.9,label=round(x_cutoff2,3)),hjust=1,color='red')+
            geom_text(aes(x=0,y=y_cutoff1*1.07,label=round(y_cutoff1,3)),hjust=0,color='red')+
            geom_text(aes(x=0,y=y_cutoff2*1.07,label=round(y_cutoff2,3)),hjust=0,color='red')+
            scale_color_manual(values=c('black','gray'))+
            xlab('Number of merges')+ylab('Distance of joining')+
            give_better_textsize_plot(8)
    #print(p2)
    
            #theme(legend.position = c(0.8, 0.2))
    p3<-ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))+
            geom_point(data=data.frame(x=x_mid, y=d_y), aes(x=x,y=y,color='derivative'))+
            geom_line(data=data.frame(x=x_mid, y=d_y), aes(x=x,y=y,color='derivative'))+
            geom_point(data=data.frame(x=x_mm, y=dd_y), aes(x=x,y=y,color='2nd_derivative'))+
            geom_line(data=data.frame(x=x_mm, y=dd_y), aes(x=x,y=y,color='2nd_derivative'))+
            geom_vline(xintercept = x_cutoff1)+geom_text(aes(x=x_cutoff1,y=max(dd_y),label=round(x_cutoff1,3)),hjust=1)+
            geom_vline(xintercept = x_cutoff2)+geom_text(aes(x=x_cutoff2,y=max(dd_y*.9),label=round(x_cutoff2,3)),hjust=1)+
            scale_color_manual(values=c('red','blue'))+
            #theme(legend.position = c(0.8, 0.2))+
            xlab('Number of merges')+ylab('derivative')+
            give_better_textsize_plot(8)
    #p3
    
    g1<-grid.arrange(p2,p3,nrow=2)

    return(list(p1=p1,
        p2=p2,
        p3=p3,g1=g1,
        distance=c(cutoff_distance1,cutoff_distance2),
        y_cutoff=c(y_cutoff1,y_cutoff2)))
    
}
