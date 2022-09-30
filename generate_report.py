import json
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import matplotlib.pyplot as plt
from inspect import signature
from xhtml2pdf import pisa
import sqlite3

def get_sample_ids(cursor):
    '''
    ascending list of sample ids

    '''
    sample_ids = []
    for row in cursor.execute(
        """
        select id 
        from analysis_report_fusions
        order by id asc
        """):
        sample_ids.append(row[0])
    
    return sample_ids


def get_table_data(cursor, processes, table_req):
    
    table_data = {}

    for process in processes:
        table_data[process] = {"header": table_req[process], "data": []}

        attributes = ",".join(table_req[process])

        for row in cursor.execute(f"select {attributes} from analysis_report_{process} order by id asc"):
            table_data[process]["data"].append(row)
    return table_data


def generate_plot(plot_data, process):
    df = pd.DataFrame(plot_data, columns=["ids", "calls"])
    df.plot(x="ids", y="calls", kind="bar")
    plt.savefig(
        f"/.mounts/labs/gsiprojects/gsi/gsiusers/jqian/Analysis_Reports/plots/{process}_plot.png",
        bbox_inches="tight"
    )

    return f"plots/{process}_plot.png"


def get_plot_data(cursor, processes):
    plots = []

    for process in processes:
        plot_data = {"ids": [], "calls" : []}
        for row in cursor.execute(f"select id, num_calls from analysis_report_{process} order by id asc"):
            plot_data["ids"].append(row[0])
            plot_data["calls"].append(row[1])
        
        plots.append(dict(process=process, plot=generate_plot(plot_data, process)))
    
    return plots

tables = {
    "mutations":[
        "id",
        "num_calls",
        "bcf_PASS_summary_records",
        "bcf_PASS_summary_SNPs",
        "bcf_PASS_summary_MNPs",
        "bcf_PASS_summary_indels",
    ], 
    "wg_structural_variants": [
        "id",
        "num_calls",
        "bcf_PASS_summary_records",
        "bcf_PASS_summary_SNPs",
        "bcf_PASS_summary_MNPs",
        "bcf_PASS_summary_indels",
    ],
    "mavis_structural_variants":[
        "id",
        "num_fusions",
    ],
    "fusions":[
        "id",
        "num_fusions",
    ],
    "copynumber_segmentation":[
        "id",
        "num_cnvs",
    ]
}

def generate_report(table_req):
    processes = ["mutations","wg_structural_variants", "fusions", "mavis_structural_variants", "copynumber_segmentation"]
    # sections = [generate_id_list, generate_table, generate_plot]

    context = {}
    
    con = sqlite3.connect("test.sl3")
    cur = con.cursor()

    context["sample_ids"] = get_sample_ids(cur)
    context["processes"] = processes
    context["table_data"] = get_table_data(cur, processes, table_req)
    # context["plots"] = get_plot_data(cur, processes)

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
    
    cur.close()
    con.close()

generate_report(tables)