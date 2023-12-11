import numpy as np
import matplotlib.pyplot as plt
import ROOT


global a
a=0.4


def epsilon_model_0(d):
     
     plasma_speed = 1
     
     height = 1    

    
     epsilon = d-(height*0.5)

     return epsilon


def epsilon_model_1(d):
     
   
     
     height = 1    

     plasma_speed = 1
     plasma_speed_constant=0.75
     
     if d<0.5:
          plasma_speed_constant = 0.75
          plasma_speed = 1
     if d>0.5:
          plasma_speed_constant = 1/0.75
          plasma_speed = 1*0.75

     epsilon = (height*0.5)*((d*(plasma_speed_constant+1))-height)/((d*((plasma_speed_constant-1)))+height)
   
     return epsilon



def epsilon_model_2(d):
     
     
     if d<0.5:
          plasma_speed = 1
          k=0.5
     if d>0.5:
          plasma_speed = 0.5
          k=1/0.5


   
     height = 1    

     epsilon = ((2*d)-height)/((2*d*(k-1))+height)*(height/2)


     return epsilon


def epsilon_model_3(d):
     
     plasma_speed = 1
     
     height = 1    

     root_1 = np.sqrt(plasma_speed**2-(4*a*d/2))
     root_2 = np.sqrt(plasma_speed**2-(4*a*(height-d)/2))
    
     epsilon = (height/2)*(-root_1+root_2)/((2*plasma_speed)-root_1-root_2)
     print(epsilon)


     return epsilon



def epsilon_model_4(d):
     
     plasma_speed = 1
     
     height = 1    

     epsilon_top = (2*d-height)*plasma_speed**2  +   2*a*(d**2-(height-d)**2)

     epsilon_bot = height*plasma_speed**2  +  2*a*(d**2+height*(height-d)**2)

     epsilon = (epsilon_top/epsilon_bot)*height*0.5


     return epsilon


def epsilon_model_5(d):
     
     plasma_speed = 1
     
     height = 1    

     D = height - d

     lim = 0.2

     t1 = calc_t_model_5(d,a,plasma_speed,lim)
     t2 = calc_t_model_5(D,a,plasma_speed,lim)

     epsilon = (height*0.5)*((t1-t2)/(t1+t2))

     return epsilon


def calc_t_model_5(d,a,v,lim):

     if d > lim or d<1-lim:
        time = d/v

     if d < lim:
        time = (d/v) + (2*a/v**3)*(((d**2)-(lim**2)))

     if d>1-lim:
        time = (d/v) + (2*a/v**3)*(((d**2)-((1-lim)**2)))

     return time

        
         

def plas_prop_time_model_0(d):
     
     plasma_speed = 1
     height = 1    

     plasma_prop_time = d/plasma_speed
  
     return plasma_prop_time

def plas_prop_time_model_1(d):
     
     
   
   
     

     if d<0.5:
          plasma_speed_constant = 0.75
          plasma_speed = 1
     if d>0.5:
          plasma_speed_constant = 1/0.75
          plasma_speed = 1*0.75
     

     height = 1    

     plasma_prop_time = (1/plasma_speed)*(d*(1-(1/plasma_speed_constant))+(height/plasma_speed_constant))
  
     return plasma_prop_time


def plas_prop_time_model_2(d):
     
     if d<0.5:
          plasma_speed = 1
          k=0.5
     if d>0.5:
          plasma_speed = 0.5
          k=1/0.5


     height = 1     

     plasma_prop_time = ((2*d/plasma_speed)*(1-(1/k)))+(height/(k*plasma_speed))

     return plasma_prop_time


def plas_prop_time_model_3(d):

          plasma_speed = 1
          
          height = 1    

          root_1 = np.sqrt(plasma_speed**2-(4*a*d/2))
          root_2 = np.sqrt(plasma_speed**2-(4*a*(height-d)/2))

          plas_prop_time = ((plasma_speed -root_1)/a )   + ((plasma_speed -root_2)/a )

          return plas_prop_time



def plas_prop_time_model_4(d):

          plasma_speed = 1
          
          height = 1   

          plas_prop_time = (d/plasma_speed) + (a/((plasma_speed)**3))*(d**2+((height - d)**2))

       
          return plas_prop_time



def plas_prop_time_model_5(d):

          plasma_speed = 1
     
          height = 1    

          D = height - d

          lim = 0.2

          t1 = calc_t_model_5(d,a,plasma_speed,lim)
          t2 = calc_t_model_5(D,a,plasma_speed,lim)

          prop_time = t1

          return prop_time







def save_plot(array_x0,array_y0,array_x1,array_y1,title,xlabel,ylabel):
    
        plt.plot(array_x0,array_y0,'.',label="linear symmetric model")
        plt.plot(array_x1,array_y1,'.',label="model 5")
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        plt.legend()
        plt.savefig("plasma_plots/"+str(title)+".png")
        plt.clf()       


def plas_prop_time_calc(model,d):


     if model ==5:
          plas_prop_time = plas_prop_time_model_5(d)
     if model ==4:
          plas_prop_time = plas_prop_time_model_4(d)
     if model ==3:
          plas_prop_time = plas_prop_time_model_3(d)
     if model ==2:
          plas_prop_time = plas_prop_time_model_2(d)
     if model ==1:
          plas_prop_time = plas_prop_time_model_1(d)
     if model ==0:
          plas_prop_time = plas_prop_time_model_0(d)

     return plas_prop_time  

def epsilon_calc(model,d):

     if model ==5:
          epsilon = epsilon_model_5(d)
     if model ==4:
          epsilon = epsilon_model_4(d)
     if model ==3:
          epsilon = epsilon_model_3(d)
     if model ==2:
          epsilon = epsilon_model_2(d)
     if model ==1:
          epsilon = epsilon_model_1(d)
     if model ==0:
          epsilon = epsilon_model_0(d)

     return epsilon         

def plot_subfunction(model):
         
          distance_array = np.linspace(0, 1, 100)
          epsilon_array = []
          plasma_prop_time_array = []
    

          for i in range(len(distance_array)):
               
               d = float(distance_array[i])
               if d>0 and d<1:
                    epsilon = epsilon_calc(model,d)
                    plasma_prop_time = plas_prop_time_calc(model,d)

                    epsilon_array.append(epsilon)
                    plasma_prop_time_array.append(plasma_prop_time)

          distance_array=  np.delete(distance_array,99)
          distance_array=  np.delete(distance_array,0)

          return epsilon_array,plasma_prop_time_array,distance_array

def plot():
          
          epsilon_array0,plasma_prop_time_array0,distance_array0 = plot_subfunction(0)
          epsilon_array1,plasma_prop_time_array1,distance_array1 = plot_subfunction(4
                                                                                    )
        
          save_plot(epsilon_array0,plasma_prop_time_array0,epsilon_array1,plasma_prop_time_array1,"Plasma Propagation Time vs. Epsilon",'eps','Plasma Propagation Time')
          save_plot(distance_array0,plasma_prop_time_array0,distance_array1,plasma_prop_time_array1,"Plasma_prop_time_vs_distance",'d1','Plasma Propagation Time')
          save_plot(epsilon_array0,distance_array0,epsilon_array1,distance_array1,"Epsilon_vs_distance",'eps','d1')


def plot2():
          t1_array = np.linspace(0, 1, 100)
          z2_array = []
          z1_array = []
          tp=1

          plasma_speed = 1
          plasma_speed_constant = 0.6
          height = 1

          for i in range(len(t1_array)):

               t1=t1_array[i]
               t2=tp-t1
               
              
               if t1>0 and t1<1:
                  
                    z2 = z_calc(t1,t2,2)
                    z2_array.append(z2)

          t1_array=  np.delete(t1_array,99)
          t1_array=  np.delete(t1_array,0)

          


          for i in range(len(t1_array)):

               t1=t1_array[i]
               t2=tp-t1
               
              
               if t1>0 and t1<1:
                  
                    z1 = z_calc(t1,t2,1)
                    z1_array.append(z1)

      

          plt.plot(t1_array,z1_array,label='z calc linear')
          plt.plot(t1_array,z2_array,label='z calc decelerating plasma')
          plt.legend()
          plt.ylabel("z")
          plt.xlabel("t1")
          plt.title("z us t model 2")
          plt.savefig("plasma_plots/"+str("z us t model 2")+".png")
          plt.clf()   





def plot_3():

     # Create a ROOT canvas
     canvas = ROOT.TCanvas("canvas", "2D Histogram", 800, 600)

     # Create a TH2F histogram
     hist = ROOT.TH2F("hist", "2D Histogram", 1000, 0, 1, 1000, 0, 1)

     hist.SetMinimum(0.5)    # Set the minimum value
     hist.SetMaximum(1) 
     # Generate some random data

     t1 = np.linspace(0, 1, 1000)
     t2 = np.linspace(0, 1, 1000)

     # Fill the histogram with data  
     for i in range(0,len(t1)):
        for j in range(0,len(t2)):
          hist.Fill(t1[i], t2[j], abs(z_calc(t1[i],t2[j], 3)))
     # print(z_calc(t1[i],t2[j], 1))

     # Set the axis titles
     hist.GetXaxis().SetTitle("t1")
     hist.GetYaxis().SetTitle("t2")
     hist.SetStats(False)



     # Draw the histogram
     hist.Draw("colz")

     # Update the canvas
     canvas.Update()

     # Save the histogram as a PNG file
     canvas.SaveAs("histogram_correction_factor.png")

     # Close the canvas
     canvas.Close()




















def correction_fector(z_cons):

    L=1
    k=0.5

    f =   (1-((k*(L/2)*(1-abs(z_cons)))))
    return f

def z_model_1(t1,t2):

     L=1
     z=0

      
     if t1>0 and t2>0:
          
           z = L/2*(t1-t2)/(t1+t2)
     
     return z
      
def z_model_2(t1,t2):

     L=1
     k=0.7
     z=0

     if t1>0 and t2>0:

        z = (L/2)*((t1-t2)/(t1+t2))*(1-((k*(L/2)*(1-abs((t1-t2)/(t1+t2))))))

     return z
      

def z_calc(t1,t2, model):    

     if model ==1:
           z = z_model_1(t1,t2)
           
     if model ==2:
           z = z_model_2(t1,t2)


     if model ==3:
          z = correction_fector(z_model_1(t1,t2))
           
     return z    











plot()

#plot2()




