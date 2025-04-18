#!/bin/bash

inputs=("sartre_bnonsat_Au_phi_ab_eAu_1.0000.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0001.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0002.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0004.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0005.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0006.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0007.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0008.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0009.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0010.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0011.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0012.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0013.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0014.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0015.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0016.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0017.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0018.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0019.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0020.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0021.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0022.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0023.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0024.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0025.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0026.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0027.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0028.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0029.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0030.eicrecon.tree.edm4eic.root"\
		"sartre_bnonsat_Au_phi_ab_eAu_1.0031.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0032.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0033.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0034.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0035.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0036.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0037.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0038.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0039.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0040.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0041.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0042.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0043.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0044.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0045.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0046.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0047.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0048.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0049.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0050.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0051.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0052.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0053.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0054.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0055.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0056.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0057.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0058.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0059.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0060.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0061.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0062.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0063.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0064.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0065.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0066.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0067.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0068.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0069.eicrecon.tree.edm4eic.root"\
        "sartre_bnonsat_Au_phi_ab_eAu_1.0070.eicrecon.tree.edm4eic.root") 
output_dir="./output"
combined_file1="${output_dir}/combined_histograms.root"
combined_file2="${output_dir}/combined_2dhistograms.root"
combined_file3="${output_dir}/combined_2dhistograms_wRES.root"
combined_file4="${output_dir}/combined_histograms_wRES.root"
combined_file5="${output_dir}/combined_histograms_wRES_cut.root"
combined_file6="${output_dir}/combined_2dhistograms_wRES_cut.root"
combined_file7="${output_dir}/combined_histograms_wCUT.root"
combined_file8="${output_dir}/combined_2dhistograms_wCUT.root"

mkdir -p "${output_dir}" # Ensure the output directory exists'

for input in "${inputs[@]}"; do
    output="${output_dir}/$(basename ${input%.*})"

    echo "Running analysis for file: ${input}"
    echo "Input edm4eic data is at [ ${input} ]"
    echo "Output at [ ${output}_output.root ]"

    root -b -q src/diffractive_vm_simple_analysis_mod.cxx+\(\"${input}\",\"${output}\"\)

    echo "Generating plots for file: ${output}_output.root"
    root -b -q macros/plot_diffractive_vm_physics_benchmark.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_diffractive_vm_resolution.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_diffractive_event_kinematics.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_2d.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_2d_wRES.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_wRES.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_wRES_cut.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_2d_wRES_cut.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_wCUT.C\(\"${output}_output.root\"\)
    root -b -q macros/plot_2d_wCUT.C\(\"${output}_output.root\"\)

    echo "Finished processing file: ${input}"
    echo
done

echo "All files processed!"

echo "Combining histograms from all processed files..."
root -b -q macros/combine_histograms.C+\(\"${output_dir}\",\"${combined_file1}\"\)
echo "Combined histograms saved to: ${combined_file1}"

echo "Generating plot for combined histograms..."
root -b -q macros/plot_diffractive_vm_combined_dsigma_dt.C
echo "Plot saved to: ./figures/combined/combined_new_method.pdf"

echo "Combining 2d histograms from all processed files..."
root -b -q macros/combine_2d.C+\(\"${output_dir}\",\"${combined_file2}\"\)
echo "Combined histograms saved to: ${combined_file2}"

echo "Generating 2d plot for combined histograms..."
root -b -q macros/plot_2d_combined.C
echo "Plot saved to: ./figures/combined/combined_2d.pdf"

echo "Combining 2d histograms with resolution from all processed files..."
root -b -q macros/combine_2d_wRES.C+\(\"${output_dir}\",\"${combined_file3}\"\)
echo "Combined 2d histograms with resolution saved to: ${combined_file3}"

echo "Generating plot for 2d combined histograms with resolution..."
root -b -q macros/plot_2d_combined_wRES.C
echo "2d Plot with resolution saved to: ./figures/combined/combined_2d_wRES.pdf"

echo "Combining histograms with resolution from all processed files..."
root -b -q macros/combine_wRES.C+\(\"${output_dir}\",\"${combined_file4}\"\)
echo "Combined histograms with resolution saved to: ${combined_file4}"

echo "Generating plot for combined histograms with resolution..."
root -b -q macros/plot_combined_wRES.C
echo "Plot with resolution saved to: ./figures/combined/combined_wRES.pdf"

echo "Combining histograms with resolution and wedge cut from all processed files..."
root -b -q macros/combine_wRES_cut.C+\(\"${output_dir}\",\"${combined_file5}\"\)
echo "Combined histograms with resolution and wedge cut saved to: ${combined_file5}"

echo "Generating plot for combined histograms with resolution and wedge cut..."
root -b -q macros/plot_combined_wRES_cut.C
echo "Plot with resolution and wedge cut saved to: ./figures/combined/combined_wRES_cut.pdf"

echo "Combining 2d histograms with resolution and wedge cut from all processed files..."
root -b -q macros/combine_2d_wRES_cut.C+\(\"${output_dir}\",\"${combined_file6}\"\)
echo "Combined 2d histograms with resolution and wedge cut saved to: ${combined_file6}"

echo "Generating 2d plot for combined histograms with resolution and wedge cut..."
root -b -q macros/plot_2d_combined_wRES_cut.C
echo "2d Plot with resolution and wedge cut saved to: ./figures/combined/combined_2d_wRES_cut.pdf"

echo "Combining histograms with wedge cut from all processed files..."
root -b -q macros/combine_wCUT.C+\(\"${output_dir}\",\"${combined_file7}\"\)
echo "Combined histograms with wedge cut saved to: ${combined_file7}"

echo "Generating plot with wedge cut for combined histograms..."
root -b -q macros/plot_combined_wCUT.C
echo "Plot with wedge cut saved to: ./figures/combined/combined_wCUT.pdf"

echo "Combining 2d histograms with wedge cut from all processed files..."
root -b -q macros/combine_2d_wCUT.C+\(\"${output_dir}\",\"${combined_file8}\"\)
echo "Combined 2d histograms with wedge cut saved to: ${combined_file8}"

echo "Generating 2d plot for combined histograms with wedge cut..."
root -b -q macros/plot_2d_combined_wCUT.C
echo "2d Plot with wedge cut saved to: ./figures/combined/combined_2d_wCUT.pdf"
