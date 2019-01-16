#helper function to return the K value at Tropicana field. Future versions of
#this function may be exported to allow for atmospheric conditiosn to be taken
#into account for Get Transverse spin.
get_K = function(){
  return(0.00515295)
}


#' Get Transverse Spin
#'
#' Get Transverse Spin takes the Trackman data that is publicly available and
#' uses it to calculate the trasverse (useful) spin. The calculations in this
#' package are based on the following paper by Dr. Allen Nathan:
#' http://baseball.physics.illinois.edu/trackman/SpinAxis.pdf
#' and all calculations were directly adapted from this workbook, also by
#' Dr. Nathan:
#' http://baseball.physics.illinois.edu/trackman/MovementSpinEfficiencyTemplate.xlsx
#'
#' @note  The current version of this function is based on atmospheric
#' conditions inside Tropicana Field in June, but future versions should
#' include paramters to specify atmospheric condiditions so as to allow for
#' more accurate spin estimates across ballparks.
#'
#' @param extension Release Extension (In Feet From Mound)
#' @param tot_spin Release Total Spin Rate (RpM)
#' @param vx0 Initial X Velocity (MpH)
#' @param vy0 Initial Y Velocity (MpH)
#' @param vz0 Initial Z Velocity (MpH)
#' @param ax X Acceleration
#' @param ay Y Acceleration
#' @param az Z Acceleration
#' @param direction Which spin component should be returned? See details
#'
#' @return If directional = TRUE, the function will return a named list with
#' total transverse spin as well as the spin in the X,Y, and Z directions,
#' otherwise it will return a scalar representing the total transverse spin.
#'
#' @details The direction argument will return the total transverse spin if
#' direction=="TRANSVERSE". "X","Y", or "Z" will return the transvers spin in
#' that direction. Direction == "ALL" will return a named list consisting of
#' the transverse spin (spinT) as well as all components (spinTx,spinTy,spinTz)
#'
#' @export
get_transverse_spin = function(extension,tot_spin,vx0,vy0,vz0
                        ,ax,ay,az,direction = "TRANSVERSE"){

  #check for validity of direction argument before doing any calculations
  if(!toupper(direction) %in% c("TRANSVERSE","X","Y","Z","ALL")){
    stop(paste0("Invalid input for direction: ",direction
                ,".\n Please select from one of: 'TRANSVERSE','X','Y','Z', "
                ,"or 'ALL'.\n For more information type ?get_transverse_spin"))
  }

  #determine release Y from extension
  yR = 60.5-extension

  #get time difference between release and y=50 ft (where v0 is measured)
  tR = (-vy0-sqrt(vy0^2-2*ay*(50-yR)))/ay

  #get component velocity at release point
  vxR = vx0+ax*tR
  vyR = vy0+ay*tR
  vzR = vz0+az*tR

  #flight time from release to y = 17/12
  tf=(-vyR-sqrt(vyR^2-2*ay*(yR-17/12)))/ay

  #calculate avg velos
  vxbar=(2*vxR+ax*tf)/2
  vybar=(2*vyR+ay*tf)/2
  vzbar=(2*vzR+az*tf)/2

  #Get total mean velocity from the vector components
  vbar =sqrt(vxbar^2+vybar^2+vzbar^2)

  #Calculate drag coefficient
  adrag = -(ax*vxbar+ay*vybar+(az+32.174)*vzbar)/vbar

  #Calculate magnus forces
  amagx = ax+adrag*vxbar/vbar
  amagy = ay+adrag*vybar/vbar
  amagz = az+adrag*vzbar/vbar+32.174
  amag = sqrt(amagx^2+amagy^2+amagz^2) # total magnus force

  #calculate lift coefficient
  CL = amag/(get_K()*vbar^2)

  #Calculate "Spin factor"
  S = 0.4*CL/(1-2.32*CL)

  #Get transverse ("useful") spin and all components
  SpinT = 78.92*S*vbar
  SpinTx = SpinT*(vybar*amagz-vzbar*amagy)/(amag*vbar)
  SpinTy = SpinT*(vzbar*amagx-vxbar*amagz)/(amag*vbar)
  SpinTz = SpinT*(vxbar*amagy-vybar*amagx)/(amag*vbar)

  #return output based on directional input
  if(toupper(direction) == "TRANSVERSE"){
    return(SpinT)
  }else if(toupper(direction) == "X"){
    return(SpinTx)
  }else if(toupper(direction) == "Y"){
    return(SpinTy)
  }else if(toupper(direction) == "Z"){
    return(SpinTz)
  }else if(direction == "ALL"){
    return(list("SpinT"=SpinT,"SpinTx"=SpinTx,"SpinTy"=SpinTy,"SpinTz"=SpinTz))

  }
}
