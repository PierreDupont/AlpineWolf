if(Sys.info()['user'] == 'pidu') { ## Pierre
  gitDir <- 'C:/My_documents/AlpineWolf'
  dataDir <- 'C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/AlpineWolf/01_Data'  
  analysisDir <- 'C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/AlpineWolf/02_Analysis'
  simDir <- 'C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/AlpineWolf/03_Simulations'
  meetDir <- 'C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/AlpineWolf/04_Meetings'
  reportDir <- 'C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/AlpineWolf/06_Report'
} else if(Sys.info()['user'] == 'virginia') {## Virginia
  gitDir <- '/Users/virginia/Dropbox/Mac/Documents/GitHub/AlpineWolf'
  dataDir <- '/Users/virginia/Dropbox/AlpineWolf/01_Data'
  analysisDir <- '/Users/virginia/Dropbox/AlpineWolf/02_Analysis'
  simDir <- '/Users/virginia/Dropbox/AlpineWolf/03_Simulations'
  meetDir <- '/Users/virginia/Dropbox/AlpineWolf/04_Meetings'
  reportDir <- '/Users/virginia/Dropbox/AlpineWolf/06_Report'
  
} else if(Sys.info()['user'] == 'richbi') {## Richard
  gitDir <- 'C:/Users/richbi/OneDrive - Norwegian University of Life Sciences/PROJECTS/Rgit/AlpineWolf'
  dataDir <- 'C:/Users/richbi/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/01_Data'  
  analysisDir <- 'C:/Users/richbi/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/02_Analysis'
  simDir <- 'C:/Users/richbi/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/03_Simulations'
  meetDir <- 'C:/Users/richbi/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/04_Meetings'
  reportDir <- 'C:/Users/richbi/Dropbox (AQEG)/AQEG Team Folder/AlpineWolf/06_Report'
  
} else stop('unknown user')

