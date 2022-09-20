import json
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import matplotlib.pyplot as plt
from inspect import signature
from xhtml2pdf import pisa
import pdfkit

def get_sample_ids(input_json):
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    sample_ids = list(data.keys())
    sample_ids.sort()

    return sample_ids

def get_table_data(input_json, table_name):
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    table_data = []

    for sample in get_sample_ids(input_json):
        table_data.append(data[sample][table_name])
        table_data[-1]["id"] = sample

    return table_data

def generate_table(input_json, table_name):

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("table.html")

    context = {
        "table_data": get_table_data(input_json, table_name),
        "process": table_name
        }

    return results_template.render(context)

    # with open("sample_table.html", "w", encoding="utf-8") as results:
    #     results.write(results_template.render(context))
    #     print("wrote to sample_table.html")

def generate_id_list(input_json):

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("sample_ids.html")

    context = {
        "sample_ids": get_sample_ids("output.json")
        }

    return results_template.render(context)

    # with open("report.html", "w", encoding="utf-8") as results:
    #     results.write(results_template.render(context))
    #     print("wrote to sample_id_.html")



def get_plot_data(input_json, plot_name):
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    plot_data = {
        "ids": [],
        "calls": []
        }
    
    plot_data["ids"] = get_sample_ids(input_json)

    for id in plot_data["ids"]:
        plot_data["calls"].append(data[id][plot_name]["num_calls"])
    
    return plot_data

#NOTE: consider generate all the plots at once??
def generate_plot(input_file, plot_name):
    df = pd.DataFrame(get_plot_data(input_file, plot_name), columns=["ids", "calls"])
    print(df)
    df.plot(x="ids", y="calls", kind="bar")
    plt.savefig(f"/home/jqian/reports/Analysis_Reports/plots/{plot_name}_plot.png", bbox_inches="tight")

    context = {
    "plot": f"plots/{plot_name}_plot.png",
    "plot_name": plot_name
    }

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("plot.html")

    return results_template.render(context)

    # with open("sample_plot.html", "w", encoding="utf-8") as results:
    #     results.write(results_template.render(context))
    #     print("wrote to sample_plot_.html")


def generate_report(input_file):
    processes = ["mutations","wg_structual_variants"]
    sections = [generate_id_list, generate_table, generate_plot]

    contents = ''

    for section in sections:
        repeat = True 
        for process in processes:
            if repeat:
                num_args = len(signature(section).parameters)
                repeat = num_args > 1
                if(num_args >1):
                    contents = contents + section(input_file, process)
                else:
                    contents = contents + section(input_file)
            else:
                break


    # contents = contents + (generate_id_list(input_file))
    # contents = contents + (generate_table(input_file, "mutations"))
    # contents = contents + (generate_plot(input_file, "mutations"))

    with open("sample_report.html", "w", encoding="utf-8") as results:
        results.write(contents)
        print("wrote to sample_report.html")

    # Convert HTML to PDF
    with open('report.pdf', "w+b") as out_pdf_file_handle:
        pisa.CreatePDF(
            src=contents,  # HTML to convert
            dest=out_pdf_file_handle)  # File handle to receive result

    # options = {
    # 'page-size': 'Letter',
    # 'margin-top': '0.35in',
    # 'margin-right': '0.75in',
    # 'margin-bottom': '0.75in',
    # 'margin-left': '0.75in',
    # 'encoding': "UTF-8",
    # 'no-outline': None,
    # 'enable-local-file-access': None
    # }
    # with open("sample_report.html") as f:
    #     pdfkit.from_file(f, "report2.pdf", options=options)

generate_report("output.json")