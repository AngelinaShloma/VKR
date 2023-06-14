from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter.messagebox import showerror
import numpy as np
import csv
import os
import pandas as pd
from method_taylor import calculate_Taylor
from rk45 import calculate_RK

filepath = "dataframe.txt"

def load_data():
    global dataframe
    dataframe = pd.read_table(filepath, sep='\s', engine='python', header=None)

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
    result.configure(state=NORMAL)
    result.delete(0.0, END)
    result.configure(state=DISABLED)

    if (os.stat(filepath).st_size == 0):
        showerror("Ошибка", "Введите данные!")
        return

    try:
        float(k1.get())
        float(k2.get())
        float(k3.get())
        float(k4.get())
        float(k5.get())
        float(k6.get())
    except TclError:
        showerror("Ошибка", "Значение параметров могут быть только числовыми!")
        return

    k = np.array([[k1.get(), k2.get(), k3.get()], [k4.get(), k5.get(), k6.get()]])
    for i in range(0, 2):
        for j in range(0, 3):
            if (k[i][j] == 0.0):
                showerror("Ошибка", "Необходимо заполнить поля ввода значения параметров!")
                return
            
    if (check_var.get() != 0 and check_var.get() != 1):
        showerror("Ошибка", "Метод не выбран!")
        return
    
    load_data()
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

def close(window_data):
    window_data.grab_release()
    window_data.destroy()
 
def click():
    window_data = Toplevel()
    window_data.title("Ввод данных")
    window_data.geometry(f"{670}x{450}+400+100")
    window_data.resizable(False, False)
    window_data.protocol("WM_DELETE_WINDOW", lambda: close(window_data))

    def add():
        time = ent_time.get()
        amount_1 = ent_amount_1.get()
        amount_2 = ent_amount_2.get()

        try:
            int(time)
            float(amount_1)
            float(amount_2)
        except ValueError:
            showerror("Ошибка", "Необходимо ввести числовые значения!")
            return
        else:
            table.insert('', 'end', values=(time, amount_1, amount_2))
            ent_time.delete(0, END)
            ent_amount_1.delete(0, END)
            ent_amount_2.delete(0, END)
            table.configure(yscrollcommand=yscrollbar.set)

    def delete():
        table.delete(*table.get_children())

    def delete_one():
        item = table.selection()[0]
        table.delete(item)

    def sort_data():
        list = [(int(table.set(k, 0)), k) for k in table.get_children("")]
        list.sort()
        for index,  (_, k) in enumerate(list):
            table.move(k, "", index)

    def save_data():
        sort_data()
        with open(filepath, "w", newline="") as f:
            csvwriter = csv.writer(f, delimiter=' ')
            for row_id in table.get_children():
                row = table.item(row_id)['values']
                csvwriter.writerow(row)
        close(window_data)
    


    frame_data = Frame(window_data)
    frame_edit = Frame(window_data)
    frame_save = Frame(window_data)
    frame_edit.pack(fill="both", padx=20, pady=10)
    frame_data.pack(fill="both", expand=1, padx=20)
    frame_save.pack(fill="both", padx=20, pady=12)

    lbl_time = Label(frame_edit, text="Момент времени:")
    lbl_time.grid(row=0, column=0, padx=5, pady=3, sticky=W)
    ent_time = Entry(frame_edit, width=15)
    ent_time.grid(row=0, column=1, padx=5, pady=3)

    lbl_amount_1 = Label(frame_edit, text="Численность 1-го вида:")
    lbl_amount_1.grid(row=1, column=0, padx=5, pady=3)
    ent_amount_1 = Entry(frame_edit, width=15)
    ent_amount_1.grid(row=1, column=1, padx=5, pady=3)

    lbl_amount_2 = Label(frame_edit, text="Численность 2-го вида:")
    lbl_amount_2.grid(row=2, column=0, padx=5, pady=3)
    ent_amount_2 = Entry(frame_edit, width=15)
    ent_amount_2.grid(row=2, column=1, padx=5, pady=3)

    btn_add = Button(frame_edit, text="Добавить", relief=GROOVE, bd=2, width=12, command=add)
    btn_add.grid(row=3, column=1, padx=5, pady=6)

    col = ('time', 'amount_1', 'amount_2')
    table = ttk.Treeview(frame_data, columns=col)
    table.column('#0', width=0)

    table.heading('#1', text='Момент времени')
    table.heading('#2', text='Численность 1-го вида')
    table.heading('#3', text='Численность 2-го вида')

    yscrollbar = Scrollbar(frame_data, orient="vertical", command=table.yview)
    table.configure(yscrollcommand=yscrollbar.set)
    table.grid(row=0, column=0)
    yscrollbar.grid(row=0, column=1, sticky=N + S)


    if os.path.isfile(filepath) == 1:
        with open(filepath) as f:
            csvread = csv.reader(f, delimiter=' ')
            for row in csvread:
                table.insert("", 'end', values=row)

    btn_save = Button(frame_save, text="Ок", relief=GROOVE, bd=2, width=10, height=1, command=save_data)
    btn_save.pack(side="right", padx=15, pady=3)

    btn_clear = Button(frame_save, text="Удалить", relief=GROOVE, bd=2, width=10, height=1, command=delete_one)
    btn_clear.pack(side="right", padx=5, pady=3)

    btn_clear_all = Button(frame_save, text="Удалить все", relief=GROOVE, bd=2, width=10, height=1, command=delete)
    btn_clear_all.pack(side="right", padx=15, pady=3)

    window_data.grab_set()

if __name__ == "__main__":
    window = Tk()
    window.title("Идентификация")
    window.geometry(f"{900}x{500}+200+80")
    window.resizable(False, False)

    titletext = Label(text="Модель Лотки-Вольтерры конкуренции двух видов", font="time 12")
    titletext.pack(side='top', pady=15)
    button_read_data = Button(text="Внести данные", command=click, relief=GROOVE, bd=2).place(x=22, y=60, width=285, height=30)
    
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

    button_start = Button(text="Идентифицировать", command=calculator, relief=GROOVE, bd=2).place(x=22, y=410, width=285, height=30)

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