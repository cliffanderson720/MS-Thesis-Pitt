# Program To Read video and Extract Frames 
import cv2
import os

  
# Function to extract frames 
def FrameCapture(path): 
      
    # Path to video file 
    vidObj = cv2.VideoCapture(path) 
  
    # Used as counter variable 
    count = 0
  
    # checks whether frames were extracted 
    success = 1
  
    while success: 
  
        # vidObj object calls read 
        # function extract frames 
        success, image = vidObj.read() 
  
        count += 1
    
        
        if (count)%10==0:
        # Saves the frames with frame-count 
            cv2.imwrite("H30\\frames\\frame%d.TIFF" % count, image) 

# Driver Code 
if __name__ == '__main__': 
  
    # Calling the function 
    FrameCapture("C:\\Users\\cdhig\\Documents\\MS-Thesis-Pitt\\Data-analysis\\Spin time\\H30\\IMG_0441.MOV") 
