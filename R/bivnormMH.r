bivnormMH<-function(rho, rho1 = 0.9, sigma = c(1.2, 1.2), steps = 1000, type = 'ind'){

    if(rho < -1 | rho > 1){
        stop("rho must be between -1 and 1")
    }

    if(steps<100)
        warning("You should really do more than 100 steps")

    targetSample.df<-data.frame(x=rep(0,steps),y=rep(0,steps))
    candidate.df<-data.frame(x=rep(0,steps),y=rep(0,steps))

    if(length(grep('^[Rr]',type))>0){
        type <- 'rw'
    }else if(length(grep('^[Ii]',type))>0){
        type <- 'ind'
    }else if(length(grep('^[Bb]',type))>0){
        type <- 'block'
    }else if(length(grep('^[Gg]',type))>0){
        type <- 'gibbs'
    }else{
        stop("Type must be one of rw, ind, block or gibbs")
    }

    x0<-c(0,0)
    x1<-c(0,0)
    mu<-c(0,0)

    if(type=='rw'){

        startValue<-c(2,2)

        sigma1<-0.5
        var1<-sigma1^2

        sigma2<-1
        var2<-sigma2^2

        k<-2*pi/1000

        u<-runif(steps)

        z1<-rnorm(steps,0,sigma1)
        z2<-rnorm(steps,0,sigma1)

        w<-1-rho^2

        targetSample.df[1,]<-startValue

        x1 <- targetSample.df[1,]
        mu <- targetSample.df[1,]
        x0 <- c(mu[1]+z1[1], mu[2]+z2[1])

        candidate.df$x <- rep(0, steps)
        candidate.df$y <- rep(0, steps)

        candidate.df[2,] <- x0

        for(n in 2:steps){
            n1 <- n-1

            x1<-targetSample.df[n1,]
            x0<-candidate.df[n,]


            canDens <- exp(-1/(2*var2*w)*(x0[1]^2-2*rho*x0[1]*x0[2] +x0[2]^2))
            curDens <- exp(-1/(2*var2*w)*(x1[1]^2-2*rho*x1[1]*x1[2] +x1[2]^2))


            if(u[n] < canDens/curDens){ ## candidate accepted
                targetSample.df[n,] <- x0
            }else{## candidate rejected
                targetSample.df[n,] <- targetSample.df[n1,]
            }

            mu <- targetSample.df[n,]
            x0 <- c(mu[1]+z1[n], mu[2]+z2[n])
            candidate.df[n+1,] <- x0

        }
    }else if(type == 'ind'){
        if(rho1 < -1 | rho > 1)
            stop("rho1 must be between -1 and 1")

        if(any(sigma<=0))
            stop("The elements of sigma must be strictly positive non-zero")

        u<-runif(steps)

        startValue <- c(2.0, 1.5)
        x1 <- startValue
        x0 <- matrix(rnorm(2*steps, rep(0,2*steps), rep(sigma,steps)),nc=2)

        x0[,2]<-rho1*x0[,1]+sqrt(1-rho1^2)*x0[,2]

        gDensInd<-function(x, rho){

            y<-x[2]
            x<-x[1]

            return(exp(-0.5/(1-rho^2)*(x^2 - 2*rho*x*y + y^2)))
        }

        qDens<-function(x, rho, rho1, sigma){
            zy<-x[2]/sigma[2]
            zx<-x[1]/sigma[1]

            return(exp(-0.5/(1-rho1^2)*(zx^2 - 2*rho*zx*zy + zy^2)))
        }

        for(n in 1:steps){
            targetSample.df[n,] <- x1
            candidate.df[n,] <- x0[n,]

            cand <- gDensInd(x0[n,], rho)
            cur <- gDensInd(x1, rho)

            qcand <- qDens(x0[n,], rho, rho1, sigma)
            qcur <- qDens(x1, rho, rho1, sigma)

            ratio <- (cand/cur)*(qcur/qcand)

            if(u[n] < ratio)
                x1 <- x0[n,]

        }
    }else if(type=='block'){
        twosteps <- 2*steps
        targetSample.df <- data.frame(x = rep(0, twosteps),
                                      y = rep(0, twosteps))
        candidate.df <- targetSample.df
        u<-runif(twosteps)

        startValue <- c(2, 1.5)
        x1<-startValue

        sx<-sqrt(1-rho^2)
        vx<-sx^2

        sigma1<-0.75
        var1<-sigma1^2

        gDensBlock<-function(x, mx, vx){
            return(exp(-0.5*(x-mx)^2/vx))
        }

        for(n in 1:steps){
            ## draw from Block 1
            mx <- rho*x1[2]

            ## draw candidate

            x0<-c(rnorm(1,mx,sigma1), x1[2])

            n1<-2*n-1

            targetSample.df[n1,]<-x1
            candidate.df[n1,]<-x0

            cand <- gDensBlock(x0[1], mx, vx)
            cur <- gDensBlock(x1[1], mx, vx)

            qcand <- gDensBlock(x0[1], mx, var1)
            qcur <- gDensBlock(x1[1], mx, var1)

            ratio<-(cand/cur)*(qcur/qcand)

            if(u[n1] < ratio)
                x1<-x0

            ## draw from block 2

            my<-rho*x1[1]

            ## draw candidate

            x0 <- c(x1[1], rnorm(1, my, sigma1))

            n2<-2*n

            targetSample.df[n2,]<-x1
            candidate.df[n2,]<-x0

            cand <- gDensBlock(x0[2], my, vx)
            cur <- gDensBlock(x1[2], my, vx)

            qcand <- gDensBlock(x0[2], my, var1)
            qcur <- gDensBlock(x1[2], my, var1)

            ratio<-(cand/cur)*(qcur/qcand)

            if(u[n2] < ratio)
                x1<-x0
        }
    }else if(type=='gibbs'){
        twosteps <- 2*steps
        targetSample.df <- data.frame(x = rep(0, twosteps),
                                      y = rep(0, twosteps))
        candidate.df <- targetSample.df
        u<-runif(twosteps)

        startValue <- c(2, 1.5)
        x1 <- startValue

        sx<-sqrt(1-rho^2)
        mx<-rho*x1[1]
        vx<-sx^2

        sigma1<-sx
        var1<-sigma1^2



        gDensGibbs <- function(x, mx, vx){
            return(exp(-0.5*(x-mx)^2/vx))
        }

        for(n in 1:steps){

            mx<-rho*x1[1]

            ## draw a candidate
            x0 <- c(rnorm(1, mx, sigma1), x1[2])

            n1<-2*n-1

            targetSample.df[n1, ] <- x1
            candidate.df[n1,] <- x0

            cand <- gDensGibbs(x0[1], mx, vx)
            cur <- gDensGibbs(x1[1], mx ,vx)

            qcand <- gDensGibbs(x0[1], mx, var1)
            qcur <- gDensGibbs(x1[1], mx ,var1)

            ratio <- (cand/cur)*(qcur/qcand)

            if(u[n1]<-ratio)
                x1 <- x0

            ## draw a candidate

            mx <- rho * x1[1]
            x0 <- c(x1[1], rnorm(1, mx, sigma1))

            n2<-2*n

            targetSample.df[n2,] <- x1
            candidate.df[n2,] <- x0

            cand <- gDensGibbs(x0[2], mx, vx)
            cur <- gDensGibbs(x1[2], mx, vx)

            qcand <- gDensGibbs(x0[2], mx, var1)
            qcur <- gDensGibbs(x1[2], mx, var1)

            ratio <- (cand/cur)*(qcur/qcand)

            if(u[n2]<ratio)
                x1 <- x0

        }
    }

    plot(y~x,data = targetSample.df,type="l")
    invisible(list(targetSample = targetSample.df))
}

