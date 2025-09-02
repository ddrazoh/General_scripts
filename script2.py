# Ali Ansaripur
from Tkinter import *
import pyperclip
import subprocess

root = Tk()
root.geometry("900x900")
pass_details = StringVar()
myList = []

def see_wifi_pass():
    global myList
    wifi_name = "Betty"  # Replace with your Wi-Fi network name
    try:
        # macOS-specific command to retrieve Wi-Fi password
        data = subprocess.check_output(['security', 'find-generic-password', '-wa', wifi_name]).decode('utf-8').strip()
        myList.append("------------------------")
        myList.append("Wifi-->" + wifi_name)
        myList.append("Password-->" + data)
        myList.append("------------------------")
    except subprocess.CalledProcessError:
        myList.append("Failed to retrieve Wi-Fi password. Make sure the Wi-Fi name is correct.")
        myList.append("------------------------")

def show_wifi_pass():
    def listToString(s):
        myStr = ""
        for ele in s:
            myStr += ele + "\n"
        return myStr
    
    myStr = listToString(myList)
    pass_details.set(myStr)

def copytoclipboard():
    password = pass_details.get()
    pyperclip.copy(password)

Label(root, text="ALI's Wifi Password Checker", font="calibri 20 bold").place(x=60, y=50)
Button(root, text=">> Initiate Process Now <<", command=see_wifi_pass).place(x=60, y=90)
Button(root, text=">> Show wifi pass details <<", command=show_wifi_pass).place(x=60, y=130)
Entry(root, textvariable=pass_details).place(width=800, height=50, x=60, y=160)
Button(root, text=">> Copy to clipboard <<", command=copytoclipboard).place(x=60, y=220)

root.mainloop()
