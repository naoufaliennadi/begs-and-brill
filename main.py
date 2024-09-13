import bbmath
import PySimpleGUI as sg
import numpy as np


# All the stuff inside your window.
layout = [[sg.Frame('',[[sg.Button("Switch mode"),sg.Text("Metric",key="unit")]])],
          [sg.Text("Pipe diameter"),sg.InputText(size=(10,1)),sg.Text("m", key="d1")],
          [sg.Text("Pipe roughness"),sg.InputText(size=(10,1)),sg.Text("m", key="d2")],
          [sg.Text("Pipe angle"),sg.InputText(size=(10,1)),sg.Text("degrees")],
          [sg.Text("")],
          [sg.Text("Oil rate"),sg.InputText(size=(10,1)),sg.Text("m^3/s", key="vr1")],
          [sg.Text("Oil density"),sg.InputText(size=(10,1)),sg.Text("kg/m^3", key="de1")],
          [sg.Text("Oil viscosity"),sg.InputText(size=(10,1)),sg.Text("Pa*s", key="v1")],
          [sg.Text("")],
          [sg.Text("Gas rate"),sg.InputText(size=(10,1)),sg.Text("m^3/s", key="vr2")],
          [sg.Text("Gas density"),sg.InputText(size=(10,1)),sg.Text("kg/m^3", key="de2")],
          [sg.Text("Gas viscosity"),sg.InputText(size=(10,1)),sg.Text("Pa*s", key="v2")],
          [sg.Text("")],
          [sg.Text("Surface tension"),sg.InputText(size=(10,1)),sg.Text("N/m", key="s")],
          [sg.Button('Ok'), sg.Button('Cancel')]]

# Create the Window
window = sg.Window('Hello Example', layout)

# starting the calculator
calc = bbmath.pgcalc()
mode = "metric"

# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()

    # if user closes window or clicks cancel
    if event == sg.WIN_CLOSED or event == 'Cancel':
        break

    if event =="Switch mode":
        calc.switch_mode()
        if mode == "metric":
            window["unit"].update("Field")
            window["d1"].update("ft")
            window["d2"].update("ft")
            window["s"].update("lbf/ft")
            window["vr1"].update("ft^3/s")
            window["vr2"].update("ft^3/s")
            window["de1"].update("lb/ft^3")
            window["de2"].update("lb/ft^3")
            window["v1"].update("cP")
            window["v2"].update("cP")
            mode = "field"
        elif mode == "field":
            window["unit"].update("Metric")
            window["d1"].update("m")
            window["d2"].update("m")
            window["s"].update("N/m")
            window["vr1"].update("m^3/s")
            window["vr2"].update("m^3/s")
            window["de1"].update("kg/m^3")
            window["de2"].update("kg/m^3")
            window["v1"].update("Pa*s")
            window["v2"].update("Pa*s")
            mode = "metric"

    if event == "Ok":
        diameter = float(values[0])
        epsilon = float(values[1])
        theta = float(values[2])
        oilcr = float(values[3])
        roho = float(values[4])
        muo = float(values[5])
        gascr = float(values[6])
        rohg = float(values[7])
        mug = float(values[8])
        sigma = float(values[9])

        oilr = oilcr / ((np.pi /4) * diameter ** 2) 
        gasr = gascr / ((np.pi /4) * diameter ** 2)
        mr = oilr + gasr
        lambdal, nfr, nvl = calc.get_dimensionless(oilr, gasr, diameter, sigma, roho)
        rohn, mun = calc.get_mix_properties(lambdal, roho, rohg, muo, mug)
        flowpattern = calc.det_flow_pattern(lambdal, nfr, theta)
        hlo = calc.get_liquid_holdup(lambdal, nfr, nvl, theta, flowpattern)
        ren = calc.get_reynold_num(lambdal, rohn, mun, diameter, mr)
        ff = calc.get_friction_factor(ren, diameter, epsilon)
        ffc = calc.correct_friction_f(lambdal, hlo, ff)
        rohs = calc.get_rohs(hlo, roho, rohg)
        pressure_gradient = calc.get_pressure_gradient(ffc, rohn, mr, diameter, rohs, theta)
        window.extend_layout(window,[[sg.Text("The pressure gradient is {:.2f} Pa/m".format(pressure_gradient))]])
window.close()
