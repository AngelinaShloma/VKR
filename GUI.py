from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
import numpy as np
import pandas as pd
from method_taylor import calculate_Taylor 
from rk45 import calculate_RK


def load_data():
    filepath = filedialog.askopenfilename()
    global dataframe
    dataframe = pd.read_table(filepath,sep='\s', engine='python', header=None)

def draw(x1, x2, x_sol1, x_sol2, t):
    ax.clear()
    ax.plot(t, x1, '-', t, x2, '-', t, x_sol1, '--', t, x_sol2, '--')
    canvas.draw()
    
def write(text):
    p = [text[0][0], text[0][1], text[0][2], text[1][0], text[1][1], text[1][2]]
    result.configure(state=NORMAL)
    result.insert(END, ', '.join(map(str, p)))
    result.configure(state=DISABLED)

def calculator():
    k = np.array([[k1.get(), k2.get(), k3.get()], [k4.get(), k5.get(), k6.get()]])
    datax = np.array([dataframe[1], dataframe[2]])
    datayears = np.array(dataframe[0])
    T = datayears.size
    x1 = np.zeros((2, T))
    x2 = np.zeros((2, T))
    if check_var.get() == 0:
        x1, x2, param = calculate_Taylor(k, datax, datayears)
        draw(datax[0], datax[1], x1, x2, datayears)
        write(param)
    if check_var.get() == 1:
        x1, x2, param = calculate_RK(k, datax, datayears)
        draw(datax[0], datax[1], x1, x2, datayears)
        write(param)

if __name__ == "__main__":
    window = Tk()
    window.title("Идентификация")
    window.geometry(f"{900}x{500}+200+80")
    window.resizable(False, False)

    titletext = Label(text="Модель Лотки-Вольтерры конкуренции двух видов", font="time 12")
    titletext.pack(side='top', pady=15)
    button_read_data = Button(text="Загрузить данные", command=load_data, relief=GROOVE, bd=2).place(x=22, y=60, width=285, height=30)

    label_k1 = Label(text='k1 :',).place(x=5, y=120, width=50, height=25)
    k1 = DoubleVar()
    entry_k1 = Entry(width=7, justify='center', textvariable=k1).place(x=55, y=120, width=50, height=25)

    label_k2 = Label(text='k2 :',).place(x=105, y=120, width=50, height=25)
    k2 = DoubleVar()
    entry_k2 = Entry(width=7, justify='center', textvariable=k2).place(x=155, y=120, width=50, height=25)
    label_k3 = Label(text='k3 :',).place(x=205, y=120, width=50, height=25)
    k3 = DoubleVar()
    entry_k3 = Entry(width=7, justify='center', textvariable=k3).place(x=255, y=120, width=50, height=25)

    label_k4 = Label(text='k4 :',).place(x=5, y=155, width=50, height=25)
    k4 = DoubleVar()
    entry_k4 = Entry(width=7, justify='center', textvariable=k4).place(x=55, y=155, width=50, height=25)
    label_k5 = Label(text='k5 :',).place(x=105, y=155, width=50, height=25)
    k5 = DoubleVar()
    entry_k5 = Entry(width=7, justify='center', textvariable=k5).place(x=155, y=155, width=50, height=25)

    label_k6 = Label(text='k6 :',).place(x=205, y=155, width=50, height=25)
    k6 = DoubleVar()
    entry_k6 = Entry(width=7, justify='center', textvariable=k6).place(x=255, y=155, width=50, height=25)


    method = Label(text="Выберите метод:").place(x=11, y=200, width=200, height=25)
       
    check_var = BooleanVar()

    method_taylor = Radiobutton(text="Метод рядов Тейлора", variable=check_var, value=0).place(x=15, y=230, width=150, height=25)
    method_rk = Radiobutton(text="Метод Рунге-Кутта", variable=check_var, value=1).place(x=7, y=260, width=150, height=25)

    button_start = Button(text="Начать идентификацию", command=calculator, relief=GROOVE, bd=2).place(x=22, y=410, width=285, height=30)

    fig = Figure(figsize=(5.4, 3.8), dpi=100)
    fig.clf()
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig)
    canvas.get_tk_widget().place(x=340, y=60)
    canvas.draw()

    method = Label(text="Результат:").place(x=340, y=450, width=60, height=25)
    result = Text(state=DISABLED)
    result.place(x=405, y=450, width=475, height=25)

    window.mainloop()