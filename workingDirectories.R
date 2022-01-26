if(Sys.info()['user'] == 'pidu') { ## Pierre
  gitDir <- 'C:/myDocuments/AlpineWolf'
  dataDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/01_Data'  
  analysisDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/02_Analysis'
  simDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/03_Simulations'
  meetDir <- 'C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/04_Meetings'
} else if(Sys.info()['user'] == 'virginia') {## Virginia
  gitDir <- '/Users/virginia/Dropbox/Mac/Documents/GitHub/AlpineWolf'
  dataDir <- '/Users/virginia/Dropbox/AlpineWolf/01_Data'
  analysisDir <- '/Users/virginia/Dropbox/AlpineWolf/02_Analysis'
  simDir <- '/Users/virginia/Dropbox/AlpineWolf/03_Simulations'
  meetDir <- '/Users/virginia/Dropbox/AlpineWolf/04_Meetings'
} else stop('unknown user')

