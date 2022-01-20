#' @title Latent state matrix (z) generation script for jags model initial values
#'
#' @description
#' \code{MakeInitsZ} returns a \code{matrix} object with the potential underlying states 
#' of the observed individuals compatible with the state-transition matrix and observation matrix in the model.
#'
#' @param y A \code{array} object containig individual observation histories (of dimensions: id*detectors*years).
#' @param STATE A \code{matrix} object representing the state-transition matrix in the model (0 or 1).
#' @param OBSERVATION A \code{matrix} object representing the observation matrix in the model (0 or 1).
#' @param f.state A \code{numeric} object representing the state of individuals when entering the data record
#' by default set to the first state in the transition matrix (1).
#' @param z.init.NA A \code{logical} object specifying whether NAs should replace the reconstructed states in the 
#' initial state matrix to avoid redundancy in the information provided to Jags (leads to an error with jags most of the time)
#' 
#' @return A \code{matrix} object with the individual state-histories z compatible with the model matrices.

MakeZ <- function( y                                ## Array of observation histories
                   , age.mx = NULL                    ## Individual age matrix (optional) 
                   , STATE                            ## State-transition matrix
                   , OBSERVATION                      ## Observation matrix
                   , f.state = 1                      ## Initial State code
                   , z.init.NA = TRUE)                ## Reconstructed states replacement by NA
{
   ### ----- 0.1.Convert 0/1 detections to 1/2 detections if needed -----
   BINARY.Y <- any(na.omit(as.vector(y))== 0)
   if(BINARY.Y){y <- y + 1}
   
   ### ----- 0.2.Convert matrix to array if needed -----
   if(is.matrix(y)){
      y.bis <- array(NA, c(dim(y)[1], 1, dim(y)[2]))
      y.bis[ ,1, ] <- y
      y <- y.bis
      }
   
   ### ----- 1.Generate a vector of first detection -----
   f <- apply(y[,1,], 1, function(x){min(which(!is.na(x)))})
   
   ### ----- 2.General attributes of the z state matrix -----
   imax <- dim(y)[1]				                      ## number of individuals
   tmax <- dim(y)[3]                                ## number of occasions	
   
   ### ----- 3.Reconstruct known states z for ids with known age -----
   z.recon <- matrix(NA, imax, tmax)	
   if(!is.null(age.mx)){z.recon <- age.mx}
   z.init <- z.recon 
   
   ### ----- 4.Create functions to read matrices -----
   events.to.states <- function(OBSERVATION, events){
      if(max(events) > dim(OBSERVATION)[2]){stop(paste("Observed event incompatible with the transition matrix for:","\n",
                                                       "    individual",i,"\n",
                                                       "    at time",t,sep=" "))}
      res <- lapply(events,function(x){which(OBSERVATION[ ,x]!=0)})
      RES <- Reduce(intersect,res)
      if(length(RES) <= 0){stop(paste("Incompatible observed events for:","\n",
                                      "    individual",i,"\n",
                                      "    at time",t,sep=" "))}
      return(RES)} # events.to.states
   
   states.to.previous.states <- function(STATE, states){
      res <- lapply(states, function(x){which(STATE[ ,x]!=0)})
      RES <- unique(unlist(res))
      return(RES)} # states.to.previous.states
   
   states.to.next.states <- function(STATE, states){
      res <- lapply(states, function(x){which(STATE[x, ]!=0)})
      RES <- unique(unlist(res))
      return(RES)} # states.to.next.states

   ### ----- 5.Sample the possible states for other individuals -----
   for(i in 1:imax){
      ## Initialize the different lists required
      st.comp.w.obs.event <- st.comp.w.next.states <- future.states <- list()
      next.states.init <- potential.z.init <- next.states.recon <- potential.z.recon <- list()
      
      ## Initialize the list of possible states (= list the states compatible with the last observation)
      future.states[[tmax]] <- events.to.states(OBSERVATION,y[i, ,tmax])
      
      ## Go through the observations backwards and list the possible states for each occasion
      ## given the current observed events and later possible states
      if(f[i] < tmax){
         for (t in (tmax-1):f[i]){
            ## List the states compatible with the observed events for each year
            st.comp.w.obs.event[[t]] <- events.to.states(OBSERVATION, y[i, ,t])
            
            ## List the states compatible with all possible future states 
            st.comp.w.next.states[[t]] <- states.to.previous.states(STATE, future.states[[t+1]])
            
            ## List the states both compatible with current observed events and possible future states
            future.states[[t]] <- intersect(st.comp.w.obs.event[[t]], st.comp.w.next.states[[t]])
            
            ## Check if the potential states list is not empty
            if(length(future.states[[t]]) <= 0){stop(paste("Incompatible potential states with future states for:","\n",
                                                "    individual",i,"\n",
                                                "    at time",t,sep=" "))}
            }#t
         }#if
      
      ## Determine potential initial states
      first.states <- intersect(f.state, future.states[[f[i]]])
      if(!is.na(z.recon[i,f[i]])){first.states <- intersect(first.states, z.recon[i,f[i]])}
      
      ## Check if the list of potential states is not empty
      if(length(first.states) <= 0){stop(paste("Proposed initial states are incompatible with initial states inferred from observations for individual",i))}
      
      ## Sample one initial state in the list of possible initial states
      z.init[i,f[i]] <- first.states[sample(length(first.states), 1)]
      z.recon[i,f[i]] <- ifelse(length(first.states) > 1, NA, z.init[i,f[i]])
      
      potential.z.recon[[f[i]]] <- first.states 
      
      for(t in (f[i]+1):tmax){
         ## List the current possible states based on the previous state (sampled or potential)
         next.states.init[[t]] <- states.to.next.states(STATE, z.init[i,t-1])
         next.states.recon[[t]] <- states.to.next.states(STATE, potential.z.recon[[t-1]])
         
         ## List the current possible states based on 1) previous state ( = next.states) and 2)future possible states ( = future states)
         potential.z.init[[t]] <- intersect(next.states.init[[t]], future.states[[t]])
         if(!is.na(z.init[i,t])){potential.z.init[[t]] <- intersect(potential.z.init[[t]], z.init[i,t])}
         potential.z.recon[[t]] <- intersect(next.states.recon[[t]], future.states[[t]])
         if(!is.na(z.recon[i,t])){potential.z.recon[[t]] <- intersect(potential.z.recon[[t]], z.recon[i,t])}
         
         ## Check if the list of potential states is not empty
         if(length(potential.z.init[[t]]) <= 0){stop(paste("Incompatible previous states with future states for:","\n",
                                                           "    individual",i,"\n",
                                                           "    at time",t,sep=" "))}
         
         ## Sample the current state in the list of possible states
         z.init[i,t] <- potential.z.init[[t]][sample(length(potential.z.init[[t]]), 1)]
         z.recon[i,t] <- ifelse(length(potential.z.recon[[t]]) > 1, NA, potential.z.recon[[t]])
         }#t
      i <- i+1
      }#i 
   
   if(z.init.NA){z.init[!is.na(z.recon)] <- NA}
   return(list("z.init.mx" = z.init, "z.reconstruct.mx" = z.recon))
} 