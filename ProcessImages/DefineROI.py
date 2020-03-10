import cv2
import numpy as np
import Image, ImageDraw
from Tkinter import *
from tkFileDialog import askopenfilename
import Image, ImageTk

def coordToMask(polygon, dims):
    img = Image.new('L', (dims[0], dims[1]), 0)
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
    return img
    
def roipoly(img2d):
    if __name__ == "__main__":
        root = Tk()
        #setting up a tkinter canvas with scrollbars
        frame = Frame(root, bd=2, relief=SUNKEN)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        xscroll = Scrollbar(frame, orient=HORIZONTAL)
        xscroll.grid(row=1, column=0, sticky=E+W)
        yscroll = Scrollbar(frame)
        yscroll.grid(row=0, column=1, sticky=N+S)
        global canvas, slice, im, x, y, polygon #for use in mouseclick function
        canvas = Canvas(frame, bd=0, xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        canvas.grid(row=0, column=0, sticky=N+S+E+W)
        xscroll.config(command=canvas.xview)
        yscroll.config(command=canvas.yview)
        frame.pack(fill=BOTH,expand=1)
        #adding the image
        slice=np.array(img2d)
        im=Image.fromarray(slice)
        rarray=np.array(im)
        #rimg=im.resize((array.shape[0]*2, array.shape[1]*2), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(im)
        canvas.create_image(0,0,image=img,anchor="nw")
        canvas.config(scrollregion=canvas.bbox(ALL))
        #mouseclick event
        x=[]
        y=[]
        polygon=[]
        canvas.bind("<Button 1>",updatepoly)
        root.mainloop()
    mask=coordToMask(polygon, slice.shape)
    del x,y,polygon,im,slice,canvas #get rid of clutter
    return mask
    
#function to be called when mouse is clicked
def updatepoly(event):
    if event.x > slice.shape[0] or event.y > slice.shape[1]: #if you click outside the image
        pass
    elif len(x) > 1:
        if abs(event.x-x[0]) < 5 and abs(event.y-y[0]) < 5: #closes polygon , is mostly aesthetic
            x.append(x[0])
            y.append(y[0])
            polygon.append(x[0])
            polygon.append(y[0])
            draw=ImageDraw.Draw(im)
            canvas.create_line(polygon, fill='red', width=1)
            canvas.update_idletasks()
        else: #normal events: updates polygon
            print (event.x,event.y)
            x.append(int(event.x))
            y.append(int(event.y))
            polygon.append(int(event.x))
            polygon.append(int(event.y))
            draw=ImageDraw.Draw(im)
            canvas.create_line([i for i in polygon], fill='red', width=1)
            canvas.update_idletasks()
    else: #prevents error message on first mouse click (line has four coordinates)
        print (event.x,event.y)
        x.append(int(event.x))
        y.append(int(event.y))
        polygon.append(int(event.x))
        polygon.append(int(event.y))
        draw=ImageDraw.Draw(im)
        canvas.create_rectangle(x[0]-1, y[0]-1, x[0]+1, y[0]+1, fill='red', width=1)
        canvas.update_idletasks()
    