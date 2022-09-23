import json
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import matplotlib.pyplot as plt
from inspect import signature
from xhtml2pdf import pisa

def get_sample_ids(data_dict):
    '''
    ascending list of sample ids

    '''
    sample_ids = list(data_dict.keys())
    sample_ids.sort()

    return sample_ids


def get_table_data(data_dict, id_list):
    
    table_data = {}

    processes = list(list(data_dict.values())[0].keys())

    for process in processes:
        table_data[process] = []

        for sample in id_list:
            table_data[process].append(data_dict[sample][process])
            table_data[process][-1]["id"] = sample

    return processes, table_data


def get_plot_data(data_dict, id_list, processes):
    plots = []

    for process in processes:
        plot_data = {"ids": id_list, "calls" : []}
        for sample_id in id_list:
            plot_data["calls"].append(data_dict[sample_id][process]["num_calls"])
        
        plots.append(dict(process=process, plot=generate_plot(id_list, plot_data, process)))
    
    return plots

def generate_plot(id_lists, plot_data, process):
    df = pd.DataFrame(plot_data, columns=["ids", "calls"])
    print(df)
    df.plot(x="ids", y="calls", kind="bar")
    plt.savefig(f"/home/jqian/reports/Analysis_Reports/plots/{process}_plot.png", bbox_inches="tight")

    return f"plots/{process}_plot.png"

    # environment = Environment(loader=FileSystemLoader("templates/"))
    # results_template = environment.get_template("plot.html")

    # return results_template.render(context)

    # with open("sample_plot.html", "w", encoding="utf-8") as results:
    #     results.write(results_template.render(context))
    #     print("wrote to sample_plot_.html")


def generate_report(input_file):
    # processes = ["mutations","wg_structual_variants"]
    # sections = [generate_id_list, generate_table, generate_plot]

    context = {}
    data = {}

    with open(input_file) as file:
        data = json.load(file)
    file.close()

    context["sample_ids"] = get_sample_ids(data)
    context["processes"], context["table_data"] = get_table_data(data, context["sample_ids"])
    context["plots"] = get_plot_data(data, context["sample_ids"], ["mutations", "wg_structual_variants"])

    with open('context.json', 'w', encoding='utf-8') as file:
        json.dump(context, file, ensure_ascii=False, indent=4)

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("base.html")

    contents = results_template.render(context)

    with open("sample_report.html", "w", encoding="utf-8") as results:
        results.write(contents)
        print("wrote to sample_report.html")

    with open('report.pdf', "w+b") as out_pdf_file_handle:
        pisa.CreatePDF(
            src=contents,
            dest=out_pdf_file_handle)

generate_report("output.json")