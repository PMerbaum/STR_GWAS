#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:09:44 2017

@author: wrheene2
"""

import argparse
import random
import math

import pandas as pd
import numpy as np
import scipy.stats as ss

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
#from itertools import chain

def get_layout():

    """General parameters for plotting"""

    layout = {"marker": ".",
              "marker_hl1": ".",
              "marker_hl2": "D",
              "marker_size": 14,
              "marker_size_hl1": 20,
              "marker_size_hl2": 20,
              "line_col": "grey",
              "line_width": 0.5,
              "font": {"fontname": "Arial", "color": "black"}}
    return layout


def colorscale(hexstr, scalefactor):
    """
    Scales a hex string by ``scalefactor``. Returns scaled hex string.

    To darken the color, use a float value between 0 and 1.
    To brighten the color, use a float value greater than 1.

    >>> colorscale("#DF3C3C", .5)
    #6F1E1E
    >>> colorscale("#52D24F", 1.6)
    #83FF7E
    >>> colorscale("#4F75D2", 1)
    #4F75D2
    """

    def clamp(val, minimum=0, maximum=255):
        if val < minimum:
            return minimum
        if val > maximum:
            return maximum
        return val

    hexstr = hexstr.strip('#')

    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr

    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)

    r = round(clamp(r * scalefactor))
    g = round(clamp(g * scalefactor))
    b = round(clamp(b * scalefactor))

    return "#%02x%02x%02x" % (r, g, b)

def get_palet(scheme):
    """Define chromosome colors
    Returns a dictionary with the colors per chromosome"""

    schemes = ["grey", "classic", "pgc", "rainbow"]
    chr = list(range(1, 23)) + ["X", "Y"]

    if scheme == "grey":
        col_even = "#AAAAAA"
        col_odd = "#CCCCCC"
        col_highlight = "#FF6600"
        col_line = "#000000"
        colors = dict(zip(chr,
                          zip([col_even, col_odd] * 12, [col_highlight] * 24)))
        colors["line"] = col_line
        return colors

    elif scheme == "classic":
        col_even = "#397DCA"
        col_odd = "#006699"
        col_highlight = "#990000"
        col_line = "#666666"
        colors = dict(zip(chr,
                          zip([col_even, col_odd] * 12, [col_highlight] * 24)))
        colors["line"] = col_line
        return colors

    elif scheme == "pgc":
        col_even = "#3458A4"
        col_odd = "#9B2C15"
        col_highlight = "#99CB99"
        col_line = "#9B2C15"
        colors = dict(zip(chr,
                          zip([col_even, col_odd] * 12, [col_highlight] * 24)))
        colors["line"] = col_line
        return colors

    elif scheme == "rainbow":
        col_chromosomes = ["#EE1100", "#F02505", "#F2390B", "#F44E11",
                           "#F66217", "#F8771D", "#FA8B23", "#FCA029",
                           "#F7AF2C", "#E4B32B", "#D0B82A", "#BDBC29",
                           "#A9C128", "#96C527", "#82CA26", "#6FCE25",
                           "#64BF32", "#5DA645", "#568D59", "#4F746C",
                           "#485B80", "#414293", "#3A29A7", "#3311BB"]
        col_highlight = ["#EE1100", "#F02505", "#F2390B", "#F44E11",
                         "#F66217", "#F8771D", "#FA8B23", "#FCA029",
                         "#F7AF2C", "#E4B32B", "#D0B82A", "#BDBC29",
                         "#A9C128", "#96C527", "#82CA26", "#6FCE25",
                         "#64BF32", "#5DA645", "#568D59", "#4F746C",
                         "#485B80", "#414293", "#3A29A7", "#3311BB"]
        col_line = "#000000"
        colors = dict(zip(chr,
                          zip(col_chromosomes, col_highlight)))
        colors["line"] = col_line
        return colors

    else:
        raise ValueError('provided scheme (%s) not in available scheme: %s'
                         % (scheme, schemes))


def get_chromosome_bounderies(build, region="autosome"):

    """Define chromosome bounderies for specified genome build.
    Returns a dataframe with length and genome position for each chromosome"""

    chromosomes = list(range(1, 23)) + ["X", "Y"]
    bounderies = pd.DataFrame({"chr": chromosomes})

    b37 = ["hg19", "37", "b37", "GRCh37"]
    b38 = ["hg38", "38", "b38", "GRCh38"]

    if region in ["aut", "autosome"]:
        region_chr = list(range(1, 23))
    elif region == "genome":
        region_chr = list(range(1, 23)) + ["X", "Y"]
    elif set(list(region)).issubset(chromosomes):
        region_chr = list(region)
    else:
        raise ValueError('region %s is not a valid genomic region: %s'
                         % (region, ["aut", "autosome", "genome"] + chromosomes))

    if str(build) in b37:
        bounderies["length"] = [249250621, 243199373, 198022430, 191154276,
                                180915260, 171115067, 159138663, 146364022,
                                141213431, 135534747, 135006516, 133851895,
                                115169878, 107349540, 102531392, 90354753,
                                81195210, 78077248, 59128983, 63025520,
                                48129895, 51304566, 155270560, 59373566]
        bounderies["end_coordinate"] = np.cumsum(bounderies["length"])
        bounderies["start_coordinate"] = [0] + list(bounderies["end_coordinate"][: -1])
        bounderies["mid_coordinate"] = (bounderies["start_coordinate"] +
                                        bounderies["end_coordinate"]) / 2
        return(bounderies[bounderies["chr"].isin(region_chr)])

    elif str(build) in b38:
        bounderies["length"] = [248956422, 242193529, 198295559, 190214555,
                                181538259, 170805979, 159345973, 145138636,
                                138394717, 133797422, 135086622, 133275309,
                                114364328, 107043718, 101991189, 90338345,
                                83257441, 80373285, 58617616, 64444167,
                                46709983, 50818468, 156040895, 57227415]
        bounderies["end_coordinate"] = np.cumsum(bounderies["length"])
        bounderies["start_coordinate"] = [0] + list(bounderies["end_coordinate"][: -1])
        bounderies["mid_coordinate"] = (bounderies["start_coordinate"] +
                                        bounderies["end_coordinate"]) / 2
        return(bounderies[bounderies["chr"].isin(region_chr)])

    else:
        raise ValueError('provided build (%s) not in available builds: %s'
                         % (build, b37 + b38))


def auto_highlight(data, threshold=None, window=250000):
    """Automatically highlight loci passing threshold"""

    print("Automatically highlight top loci...")

    if not threshold:
        threshold = 0.05 / data.shape[0]

    hl1_list = list()
    hl2_list = list()
    n_gws = data.loc[data.p < threshold].shape[0]

    while n_gws > 0:
        # highlight top hit
        i = data.loc[data.p == data.p.min()].index[0]
        hl2_list.extend(data.loc[i,"label"])
        # highlight SNPs within locus
        tophit_chr = int(data.loc[i,"chr"])
        tophit_bp = int(data.loc[i,"bp"])
        i = data.loc[(data.chr == tophit_chr) &
                     (data.bp > tophit_bp - window) &
                     (data.bp < tophit_bp + window)].index
        hl1_list.extend(data.loc[i,"label"])
        # remove highlighted SNPs
        data = data.drop(i)
        n_gws = data.loc[data.p < threshold].shape[0]

    return(hl1_list, hl2_list)


def prepare_data(data, build, region, scheme, highlight1=list(),
                 highlight2=list(), autohighlight=False, threshold=None):

    """Prepare dataframe that can be used to plot"""

    # get chromosome bounderies and colors
    palet = get_palet(scheme)
    bounderies = get_chromosome_bounderies(build, region)

    # set threshold
    if not threshold:
        threshold = 0.05 / data.shape[0]

    # set coordinates
    chr_to_start = dict(zip(bounderies.chr, bounderies.start_coordinate))
    data["coordinate"] = np.add([chr_to_start[x] for x in data["chr"]],
                                data["bp"])

    # automatically highlight loci passing threshold
    if autohighlight:
       (highlight1,highlight2) = auto_highlight(data=data, threshold=threshold)

    # set highlight 1
    print("Set highlight 1")
    data["hl1"] = np.in1d(np.array(data.label), np.array(highlight1))
    # i = data.loc[data.p < threshold].index
    # data.loc[i, "h1"] = True

    # set highlight 2 (only apply to top hit within locus with identical name)
    print("Set highlight 2")
    data["hl2"] = False
    for x in highlight2:
        locus = data.loc[data.label == x]
        i = locus.loc[locus.p == np.min(locus.p)].index
        data.loc[i, "hl2"] = True

    # set color
    print("Set colors")
    data["color"] = None
    for chr in bounderies.chr:
        i = data.loc[data.chr == chr].index
        h1 = data.loc[(data.chr == chr) & (data.hl1)].index
        h2 = data.loc[(data.chr == chr) & (data.hl2)].index
        data.loc[i, "color"] = palet[chr][0]
        data.loc[h1, "color"] = palet[chr][1]
        data.loc[h2, "color"] = palet[chr][1]

    # set edgecolor
    print("Set edge colors")
    data["edgecolor"] = data["color"][:]
    for chr in bounderies.chr:
        i = data.loc[(data.chr == chr) & (data.hl2)].index
        data.loc[i, "edgecolor"] = colorscale(palet[chr][1], 0.7)

    # limit the p-values
    print("Set upper limits of p-value")
    if min(data["p"]) < 1e-50:
        print("WARNING: extremely small p-values found, set to 1e-50")
        data["p"] = pd.Series([x if x > 1e-50 else 1e-50 for x in data["p"]])

    return data


def miami_plot(data1, data2, label1, label2, build, region,
                          filename, format="png", dpi=150, size=(15, 7.5),
                          threshold=None):

    """Make a miami plot"""

    layout = get_layout()
    bounderies = get_chromosome_bounderies(build, region)

    if not threshold:
        bar = - np.log10(0.05 / data1.shape[0])
    else:
        bar = - np.log10(threshold)

    outfile = str(filename + "_miamiplot." + format)
    ylim = np.max(np.append(np.append(- np.log10(data1["p"]), bar), -np.log10(data2["p"])))
    #ylim = 50

    fig, ax = plt.subplots(figsize=size)
    # set the axis and lines
    ax.set_ylabel("$-\log_{10}(P)$", size=10)
    ax.set_xlabel("chromosome", size=10, fontdict=layout["font"])
    ax.set_ylim([1.1 * -ylim, 1.1 * ylim])
    ax.set_ylim([-ylim, ylim])
    ax.set_yticks(list(range(-math.ceil(ylim), math.ceil(ylim)+1, 2)))
    ax.set_yticklabels(np.abs(list(range(-math.ceil(ylim), math.ceil(ylim)+1, 2))))
    ax.axes.get_xaxis().set_ticks([])
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.plot((0, np.max(data1["coordinate"])), (bar, bar), 'k--',
            c=layout["line_col"], linewidth=layout["line_width"])
    ax.plot((0, np.max(data1["coordinate"])), (-bar, -bar), 'k--',
            c=layout["line_col"], linewidth=layout["line_width"])
    for i in range(bounderies.shape[0]):
        ax.text(bounderies.mid_coordinate[i], 0, bounderies.chr[i],
                    ha="center", va="center", fontdict=layout["font"])

    # plot the data
    ax.scatter(data1["coordinate"], -np.log10(data1["p"]),
               c=data1["color"], s=layout["marker_size"], marker=layout["marker"])
    ax.scatter(data2["coordinate"], np.log10(data2["p"]),
               c=data2["color"], s=layout["marker_size"], marker=layout["marker"])
    # plot the loci to highlight
    data1_hl1 = data1.loc[data1.hl1]
    ax.scatter(data1_hl1["coordinate"], -np.log10(data1_hl1["p"]),
               c=data1_hl1["color"], s=layout["marker_size_hl1"],
               marker=layout["marker_hl1"])
    data2_hl1 = data2.loc[data2.hl1]
    ax.scatter(data2_hl1["coordinate"], np.log10(data2_hl1["p"]),
               c=data2_hl1["color"], s=layout["marker_size_hl1"],
               marker=layout["marker_hl1"])
    # plot the top hits to highlight
    data1_hl2 = data1.loc[data1.hl2]
    ax.scatter(data1_hl2["coordinate"], -np.log10(data1_hl2["p"]),
               c=data1_hl2["color"], s=layout["marker_size_hl2"],
               marker=layout["marker_hl2"], edgecolor=data1_hl2["edgecolor"])
    data2_hl2 = data2.loc[data2.hl2]
    ax.scatter(data2_hl2["coordinate"], np.log10(data2_hl2["p"]),
               c=data2_hl2["color"], s=layout["marker_size_hl2"],
               marker=layout["marker_hl2"], edgecolor=data2_hl2["edgecolor"])
    # add the labels to the top hits
    for i in range(data1_hl2.shape[0]):
        ax.text(np.array(data1_hl2.coordinate)[i], ylim * 0.025 + -np.log10(np.array(data1_hl2.p)[i]),
                np.array(data1_hl2.label)[i], ha="center", fontdict=layout["font"], fontsize=8)
    # add labels to the datasets
    ax.text(0, ylim * 0.9, label1, fontdict=layout["font"], va="top", ha="left")
    ax.text(0, - ylim * 0.9, label2, fontdict=layout["font"], va="bottom", ha="left")
    # save the file
    fig.savefig(outfile, dpi=dpi)


def manhattan_plot(data,  build, region, filename, main_title="Manhattan plot",
                   format="png", dpi=150, size=(15, 7.5), threshold=None):

    """Make a regular manhattan plot"""

    if not threshold:
        bar = - np.log10(0.05 / data.shape[0])
    else:
        bar = - np.log10(threshold)

    layout = get_layout()
    bounderies = get_chromosome_bounderies(build, region)

    outfile = str(filename + "_manhattanplot." + format)
    ylim = np.max(np.append(list(- np.log10(data["p"])), bar))
    # ylim = 50
    data_hl1 = data.loc[data.hl1]
    data_hl2 = data.loc[data.hl2]

    fig, ax = plt.subplots(figsize=size)
    # set the axis and lines
    ax.set_ylim([0, 1.1 * ylim])
    #ax.set_ylim([0, ylim])
    # plot the data
    if format == "png":
        ax.axis('off')
        ax.scatter(data["coordinate"], -np.log10(data["p"]),
                   c=data["color"], s=layout["marker_size"], marker=layout["marker"])
        # plot the loci to highlight
        ax.scatter(data_hl1["coordinate"], -np.log10(data_hl1["p"]),
                   c=data_hl1["color"], s=layout["marker_size_hl1"],
                   marker=layout["marker_hl1"])
        # plot the top hits to highlight
        ax.scatter(data_hl2["coordinate"], -np.log10(data_hl2["p"]),
                   c=data_hl2["color"], s=layout["marker_size_hl2"],
                   marker=layout["marker_hl2"], edgecolor=data_hl2["edgecolor"])

    # add the labels to the top hits
    if format == "pdf":
        ax.set_ylabel("$-\log_{10}(P)$", size=10)
        ax.set_xlabel("chromosome", size=10, fontdict=layout["font"])
        ax.set_xticks(bounderies.mid_coordinate)
        ax.set_xticklabels(bounderies.chr)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.plot((0, np.max(data["coordinate"])), (bar, bar), 'k--',
                 c=layout["line_col"], linewidth=layout["line_width"])
        ax.set_title(main_title, fontdict=layout["font"], fontsize=15)
        for i in range(data_hl2.shape[0]):
            ax.text(np.array(data_hl2.coordinate)[i], ylim * 0.02 + -np.log10(np.array(data_hl2.p)[i]),
                    np.array(data_hl2.label)[i], ha="center", fontdict=layout["font"], fontsize=8)
        # plot the loci to highlight
        ax.scatter(data["coordinate"], -np.log10(data["p"]), # added 
                   c=data["color"], s=layout["marker_size"], marker=layout["marker"]) #added 
        ax.scatter(data_hl1["coordinate"], -np.log10(data_hl1["p"]),
                   c=data_hl1["color"], s=layout["marker_size_hl1"],
                   marker=layout["marker_hl1"])
        # plot the top hits to highlight
        ax.scatter(data_hl2["coordinate"], -np.log10(data_hl2["p"]),
                   c=data_hl2["color"], s=layout["marker_size_hl2"],
                   marker=layout["marker_hl2"], edgecolor=data_hl2["edgecolor"])

    # save the file
    fig.savefig(outfile, dpi=dpi)


def qq_plot(p, ncases, ncontrols, filename, main_title="QQ plot",
            format="png", dpi=150):

    """Make a QQ plot"""

    # threshold the extremely small p-values
    if min(p) < 1e-10:
        print("WARNING: small p-values found, set to 1e-10")
        p = pd.Series([x if x > 1e-10 else 1e-10 for x in p])

    # define the parameters
    print("calculate observed and expected p-values")
    obs = sorted(-np.log10(p), reverse=True)
    n_obs = len(obs)
    exp = -np.log10([x/n_obs for x in range(1,n_obs+1)])
    print("observed " + str(n_obs) + " p-values")

    # define the plotting parameters
    colr_ci	= "#7491c1"
    colr_null = "#003366"
    pointsize = 15
    pointmarker = "."

    # calculate 95% confidence intervals
    print("calculate 95% confidence interval")
    a = list(range(1, n_obs + 1))
    b = [n_obs +1 - x for x in a]
    c95 = ss.beta.ppf(0.95, a, b)
    c05 = ss.beta.ppf(0.05, a, b)

    # make ci_polygon
    print("draw 95% confidence interval")
    c05_polygon = list(zip(exp, -np.log10(c05)))
    c95_polygon = list(zip(exp, -np.log10(c95)))
    ci_polygon = Polygon(c05_polygon + list(reversed(c95_polygon)))

    # calculate QQ values
    print("transform p-values in z-values")
    stats = 0.5 * p
    z = ss.norm.ppf(1 - stats)
    print("calculate lambda:")
    lambda_raw = np.median(z ** 2) / ss.chi2.ppf(0.5, df=1)
    print("\t" + str(lambda_raw))
    lambda_gc = np.around(lambda_raw, 3)
    if np.sum([ncases, ncontrols]) > 0:
        print("calculate lambda1000:")
        lambda_1000 = np.around(1 +
                                (lambda_raw - 1) * (1/ncases + 1/ncontrols) /
                                (1/1000 + 1/1000), 3)
        print("\t" + str(lambda_1000))
    else:
         lambda_1000 = "NA"

    print("draw QQ-plot")
    outfile = str(filename + "_qqplot." + format)
    xmax = np.max(exp) * 1.02
    ymax = np.max(obs) * 1.02
    ymax = np.max(np.concatenate([obs, -np.log10(c05)])) * 1.02
    fig, ax = plt.subplots()
    # set the axis and lines
    ax.set_xlabel("Expected " + "$-\log_{10}(P)$", size=10)
    ax.set_ylabel("Observed " + "$-\log_{10}(P)$", size=10)
    ax.set_xlim([0, xmax])
    ax.set_ylim([0, ymax])
    if format == "pdf":
        ax.plot((0, np.max(exp)), (0, np.max(exp)), ls="--", color=colr_ci, alpha=0.5, linewidth=1) #added
        ax.scatter(exp, obs, c=colr_null, s=pointsize, marker=pointmarker) #added
        # add the text
        legend_text = "$\lambda_{gc} =$" + " " + str(lambda_gc) + "\n" + "$\lambda_{1000} =$" + " " + str(lambda_1000)
        ax.text(0.2, 0.85*ymax, legend_text)
        # add the title
        ax.set_title(main_title, fontsize=15)

    # plot the data
    if format == "png":
        # plot the CI polygon
        #ci_polygon_patch = PolygonPatch(ci_polygon, color=colr_ci, alpha=0.2)
        #ax.add_patch(ci_polygon_patch)
        legend_text = "$\lambda_{gc} =$" + " " + str(lambda_gc) + "\n" + "$\lambda_{1000} =$" + " " + str(lambda_1000) #added
        ax.text(0.2, 0.85*ymax, legend_text) #added
        # plot the diagnoal
        ax.plot((0, np.max(exp)), (0, np.max(exp)), ls="--", color=colr_ci, alpha=0.5, linewidth=1)
        ax.scatter(exp, obs, c=colr_null, s=pointsize, marker=pointmarker)

    # save the plot
    fig.savefig(outfile, dpi=dpi)


def simulate_data(n=30000, scheme="grey", threshold=5e-7):

    """Simulate summary statistics"""

    bounderies = get_chromosome_bounderies("b37", "autosome")
    palet = get_palet(scheme)

    # simulate genome coordinates and find chromosomes
    coordinates = np.sort(random.sample(
                          range(np.max(bounderies["end_coordinate"])), n))
    chrs = np.digitize(coordinates, [0] + list(bounderies["end_coordinate"]))
    labels = ["rs" + str(s) for s in coordinates]

    # simulate p-values with some inflation
    betas = ss.norm(0, 1.5).rvs(n)
    null = ss.norm(0, 1)
    pvalues = [2 * null.cdf(- abs(x)) for x in betas]

    # format  dataframe
    data = pd.DataFrame({"chr": chrs, "coordinate": coordinates,
                         "label": labels, "p": pvalues})
    data["hl1"] = np.where(data["p"] < threshold, 1, 0)
    data["color"] = [palet[data.chr[x]][1] if data.hl1[x] == 1 else
                     palet[data.chr[x]][0] for x in range(n)]
    data["size"] = np.where(data["hl1"] == 1, 4, 2)

    return data


def main():

    parser = argparse.ArgumentParser(description="Plot summary statistics: manhattan plot, miami plot to compare two traits or QQ plot")

    # define the input group:
    group = parser.add_argument_group("input files")
    group.add_argument("--infile", type=str, action="store", required=True,
                       help="File with summary statistics with field chr, bp, p, label")
    group.add_argument("--infile2", type=str, action="store",
                        default=None,
                        help="File with summary statistics with field chr,bp, p, label for second dataset (for mirrored manhattan plot)")
    group.add_argument("--highlight1", type=str, action="store",
                        help="File with list of labels to highlight with different color and size, use for loci")
    group.add_argument("--highlight2", type=str, action="store",
                        help="File with list of labels to highlight with different symbol and textlabel, use for top-hit")

    group = parser.add_argument_group("file description")
    group.add_argument("--label1", type=str, action="store",
                        help="name of dataset 1")
    group.add_argument("--label2", type=str, action="store",
                        help="name of dataset 2 only for miami plot")
    group.add_argument("--ncases", type=int, action="store",
                        default=0, help="Number of cases to calculate lambda10000 in QQplot")
    group.add_argument("--ncontrols", type=int, action="store",
                        default=0, help="Number of controls to calculate lambda10000 in QQplot")
    group.add_argument("--chr-col", type=str, dest="chr_col",
                       default="chr", help="Column name for chromosome, default = chr")
    group.add_argument("--bp-col", type=str, dest="bp_col",
                       default="bp", help="Column name for basepair, default = bp")
    group.add_argument("--p-col", type=str, dest="p_col",
                       default="p", help="Column name for p-value, default = p")
    group.add_argument("--label-col", type=str, dest="label_col",
                       default="label", help="Column name for label, default = label")
    group.add_argument("--sep", type=str, action="store", default="\t",
                        help="field separator")

    group = parser.add_argument_group("plotting options")
    group.add_argument("--type", type=str, action="store",
                        choices=["manhattan", "miami", "qq"], help="type of plot")
    group.add_argument("--autohighlight", type=bool, action="store",
                       default=False,
                       help="Automatically highlight loci passing threshold")
    group.add_argument("--threshold", type=float, action="store",
                        help="P-value threshold for horizontal line in manhattan plot.")
    group.add_argument("--build", type=str, action="store",
                        default="b37", help="Build of the human genome to calculate chromosome lengths, default = b37")
    group.add_argument("--region", type=str, action="store",
                        default="autosome", help="region of the genome to plot, default = autosome, alterative = genome that includes X and Y (not fully tested)")
    group.add_argument("--colors", metavar="STRING", type=str, action="store",
                        choices=["grey", "classic", "pgc", "rainbow"],
                        default="classic", help="Color scheme, default = classic (blue with red highlights).")

    group = parser.add_argument_group("output files")
    group.add_argument("--out", metavar="FILE", type=str, action="store", required=True,
                        help="Prefix of the ouput file.")

    group = parser.add_argument_group("output file options")
    group.add_argument("--format", metavar="format", type=str, action="store",
                        default="png", help="Figure file format default = png.")
    group.add_argument("--size", metavar="size", type=tuple, action="store",
                        default=(15, 7.5), help="Figure file size default = 15 x 7.5")
    group.add_argument("--dpi", metavar="dpi", type=int, action="store",
                        default=150, help="Figure resolution, default = 150.")


    args = parser.parse_args()

    # print arguments:
    print("\n")
    print('{0:15s} {1:10s}'.format("Argument", "Value"))
    print('{0:15s} {1:10s}'.format("--------", "-----"))
    for arg, value in sorted(vars(args).items()):
        print('{0:15s} {1:10s}'.format(arg, str(value)))
    print("\n\n")

    # recode arguments as capital variables, easy to recognize
    INFILE1 = args.infile
    INFILE2 = args.infile2 if args.infile2 else None
    LABEL1 = args.label1
    LABEL2 = args.label2
    HIGHLIGHT1 = args.highlight1 if args.highlight1 else None
    HIGHLIGHT2 = args.highlight2 if args.highlight2 else None
    TYPE = args.type
    AUTOHIGHLIGHT = args.autohighlight
    NCASES = args.ncases
    NCONTROLS = args.ncontrols
    OUT = args.out
    THRESHOLD = args.threshold if args.threshold else None
    BUILD = args.build
    REGION = args.region
    COLORS = args.colors
    FORMAT = args.format
    SIZE = args.size
    DPI = args.dpi

    print("Reading datafile " + INFILE1 + " (" + LABEL1 + ")...")
    if INFILE1[-3:] == ".gz":
        data1 = pd.read_csv(INFILE1, sep=args.sep, compression="gzip", na_values="NA")
    else:
        print("Consider compressing summary statistics")
        data1 = pd.read_csv(INFILE1, sep=args.sep, na_values="NA")
    data1.rename(columns={args.chr_col: "chr", args.bp_col: "bp",
                          args.p_col: "p", args.label_col:"label"}, inplace=True)
    data1 = data1.loc[pd.notnull(data1["p"]),]

    if INFILE2:
        print("Reading datafile " + INFILE2 + " (" + LABEL2 + ")...")
        if INFILE2[-3:] == ".gz":
            data2 = pd.read_csv(INFILE2, sep=args.sep, compression="gzip", na_values="NA")
        else:
            print("Consider compressing summary statistics")
            data2 = pd.read_csv(INFILE2, sep=args.sep, na_values="NA")
        data2.rename(columns={args.chr_col: "chr", args.bp_col: "bp",
                              args.p_col: "p", args.label_col:"label"}, inplace=True)
        data2 = data2.loc[pd.notnull(data2["p"]),]


    if not HIGHLIGHT1 and not AUTOHIGHLIGHT:
        print("No labels provided to highlight loci")
        highlight1 = list()
    elif not AUTOHIGHLIGHT:
        with open(HIGHLIGHT1) as f:
            highlight1 = f.read().splitlines()
    else:
        highlight1 = list()

    if not HIGHLIGHT2 and not AUTOHIGHLIGHT:
        print("No labels provided to highlight top hits")
        highlight2 = list()
    elif not AUTOHIGHLIGHT:
        with open(HIGHLIGHT2) as f:
            highlight2 = f.read().splitlines()
    else:
        highlight2 = list()

    if TYPE == "manhattan":
        print("munge data...")
        data_ready = prepare_data(data=data1, build=BUILD, region=REGION,
                                  scheme=COLORS, highlight1=highlight1,
                                  highlight2=highlight2,
                                  autohighlight=AUTOHIGHLIGHT,
                                  threshold=THRESHOLD)
        print("plot manhattan plot to " + OUT + "_manhattanplot." + FORMAT)
        manhattan_plot(data=data_ready, build=BUILD, region=REGION,
                       filename=OUT, main_title=LABEL1, format=FORMAT, dpi=DPI,
                       size=SIZE, threshold=THRESHOLD)

    if TYPE == "miami":
        print("munge dataset 1...")
        data1_ready = prepare_data(data=data1, build=BUILD, region=REGION,
                                   scheme=COLORS, highlight1=highlight1,
                                   highlight2=highlight2,
                                   autohighlight=AUTOHIGHLIGHT,
                                   threshold=THRESHOLD)
        print("munge dataset 2...")
        data2_ready = prepare_data(data=data2, build=BUILD, region=REGION,
                                   scheme=COLORS, highlight1=highlight1,
                                   highlight2=highlight2,
                                   autohighlight=AUTOHIGHLIGHT,
                                   threshold=THRESHOLD)
        print("plot miami mirror plot to " + OUT + "_miamiplot." + FORMAT)
        miami_plot(data1=data1_ready, data2=data2_ready,
                              label1=LABEL1, label2=LABEL2,
                              build=BUILD, region=REGION, filename=OUT,
                              format=FORMAT, dpi=DPI, size=SIZE,
                              threshold=THRESHOLD)

    if TYPE == "qq":
        print("plot QQ-plot to " + OUT + "_qqplot." + FORMAT)
        qq_plot(p=data1["p"], ncases=NCASES, ncontrols=NCONTROLS, main_title=LABEL1,
                filename=OUT, format=FORMAT, dpi=DPI)


if __name__ == "__main__":
    main()
