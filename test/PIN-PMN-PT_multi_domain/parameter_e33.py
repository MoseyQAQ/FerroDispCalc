from Main import AVG_disp_polar, Output

"""此代码为刘仕老师课题组私有，请不要将其分享到GitHub上，谢谢"""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""初始输入模块"""
input_path = "./traj.lammpstrj"
Number_atoms = 16*16*16*5
Number_cell = [16, 16, 16]
timestep_selected = 1
criteria = 4
switch = 1      #switch=1代表打开强制配位数锁定，switch=0代表采用自适应配位数
charge_A = {2: 3.502}
charge_B = {9: 5.434, 11:5.434, 12:5.434, 13:5.434}
charge_X = -8.936/3
"""输出文件地址"""
output_path_AVGcell_coord = "./avg_cell_coord"
output_path_AVGdisp = "./avg_disp"
output_path_id = "./id"
output_path_AVGpolar_per_cell = "./polar_per_cell"
output_path_AVGpolar_AandB = "./polar_AandB"

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""函数执行模块"""
function = AVG_disp_polar(input_path=input_path, Number_atoms=Number_atoms, Number_cell=Number_cell,
                          timestep_selected=timestep_selected, criteria=criteria,
                          charge_A=charge_A, charge_B=charge_B, charge_X=charge_X, switch=switch)
cell_avg, coord_avg = function.ALLcoord_avg()
disp_fin, corrd_number_fin, record_avg_fin, id_AX, id_BX, error = function.get_atomsdisplacment(cell_avg, coord_avg)
Polar_cell = function.get_totPolar_cell(cell_avg, coord_avg)
PolarA, PolarB, avg_PolarA, avg_PolarB, avg_AllPolar = function.Polar_eachAtoms(cell_avg, disp_fin)
"""文件输出模块"""
output = Output(Number_atoms=Number_atoms, Number_cell=[16, 16, 16])
#output.output_avg_cell_coord(output_path_AVGcell_coord, cell_avg, coord_avg)
output.output_disp(output_path_AVGdisp, disp_fin, corrd_number_fin, record_avg_fin, error)
#output.output_id(output_path_id, id_AX, id_BX)
output.output_polar_cell(output_path_AVGpolar_per_cell, Polar_cell)
#output.output_Polar_eachAtoms(output_path_AVGpolar_AandB, PolarA, PolarB, avg_PolarA, avg_PolarB, avg_AllPolar)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
