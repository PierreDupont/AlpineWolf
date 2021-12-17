if(Sys.info()['user'] == 'pidu') { ## Pierre
  gitDir <- 'C:/myDocuments/AlpineWolf'
  dataDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/01_Data'  
  analysisDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/02_Analysis'
  simDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/03_Simulation'
} else if(Sys.info()['user'] == 'virginia') {## Virginia
  gitDir <- '.../AlpineWolf'
  dataDir <- '.../Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/01_Data'  
  analysisDir <- '.../Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/02_Analysis'
  simDir <- '.../Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/03_Simulation'
} else stop('unknown user')

