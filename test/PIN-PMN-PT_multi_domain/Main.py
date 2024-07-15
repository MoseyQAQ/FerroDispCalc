import numpy as np

class __Core__:
    """核心函数，此代码为刘仕老师课题组私有，请不要将其分享到GitHub上，谢谢"""
    def __init__(self, input_path, Number_atoms, Number_cell, timestep_selected, criteria, charge_A, charge_B, charge_X, switch):
        self.input_path = input_path
        self.Number_atoms = Number_atoms
        self.timestep_selected = timestep_selected
        self.criteria = criteria * 0.85
        self.charge_A = charge_A
        self.charge_B = charge_B
        self.charge_X = charge_X
        self.switch = switch    #switch=1代表打开强制配位数锁定，switch=0代表采用自适应配位数

        self.nx = Number_cell[0]  # Number of units in x direction
        self.ny = Number_cell[1]  # Number of units in y direction
        self.nz = Number_cell[2]  # Number of units in z direction
        self.cellAll = Number_cell[0] * Number_cell[1] * Number_cell[2]

        self.N_A = self.nx * self.ny * self.nz  # Number of A site atoms
        self.N_B = self.N_A  # Number of B site atoms
        self.N_X = 3 * self.N_A  # Number of X site atoms

        self.ff = 1.602176E-19 * 1.0E-10 * 1.0E30  # This converts polarization to C/m^2
        self.dim = 3    #数据的维度，默认是3维的
        self.error = None

        self.num_files, self.type_index, self.coord_ini, self.id = self.__para_int()     # 初始化内部参数

    def __para_int(self):
        """
            这个是内部参数初始化函数
            :return: 总时间步数num_files，原子类型type_index，初始文件的原子坐标coord_ini
        """
        str_lens = 9  # traj.lammpstrj文件开头字符串长度，也就是从"ITEM: TIMESTEP"到"ITEM: ATOMS id type x y z"的字符串长度

        with open(self.input_path, 'r') as f:
            infile = f.readlines()  # 总行数
        num_files = len(infile) // (self.Number_atoms + str_lens)  # 记录的不同时间步数目：总行数//每一时间步TIMESTEP占用的行数
        f.close()
        # num_files代表运行的总步数

        with open(self.input_path, 'r') as f:
            "寻找最初的原子坐标"
            skip_step = 9  # 跳过前9行文字描述，直接转到初始坐标处
            for j in range(skip_step):
                f.readline()
            type_index, coord_ini, id_index = self.__readatoms(f, self.Number_atoms)
        f.close()
        return num_files, type_index, coord_ini, id_index

    def __readatoms(self, f, natoms):
        """
        读取原子的坐标和对应的索引的初始化函数，如果lammps原子坐标输出文件形式
        不是按照[id type x y z]形式输出的，则需要自行修改这个函数
        :param f: 输入文件
        :param natoms: 体系原子总数
        :return: 返回原子坐标矩阵和对应的索引
        """
        type_index = np.zeros(natoms)
        id_index = np.zeros(natoms)
        coord = np.zeros((natoms, 3))
        for i in range(natoms):
            line = f.readline().split()  # line的形式[id type x y z]为代表的数据
            id_index[i] = int(line[0])
            type_index[i] = float(line[1])  # line[1]代表type
            tmp = [float(x) for x in line[2:]]  # line[2:]代表[x y z]
            coord[i, 0] = tmp[0]
            coord[i, 1] = tmp[1]
            coord[i, 2] = tmp[2]
        return type_index, coord, id_index

    def __readcell(self, f):
        """
        读取晶格常数的初始化函数
        :param f: 输入文件
        :return: 晶格常数矩阵
        """
        line = f.readline().split()  # 读取x坐标
        line = [float(x) for x in line]
        xlo_bound = line[0]
        xhi_bound = line[1]
        xy = line[2]
        line = f.readline().split()  # 读取y坐标
        line = [float(x) for x in line]
        ylo_bound = line[0]
        yhi_bound = line[1]
        xz = line[2]
        line = f.readline().split()  # 读取z坐标
        line = [float(x) for x in line]
        zlo_bound = line[0]
        zhi_bound = line[1]
        yz = line[2]
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound
        xx = xhi - xlo
        yy = yhi - ylo
        zz = zhi - zlo
        cell = np.zeros((3, 3))
        cell[0, :] = [xx, 0, 0]
        cell[1, :] = [xy, yy, 0]
        cell[2, :] = [xz, yz, zz]
        return cell

    def __data_wash(self, coord_tmp, coord_ini, cell_tmp):
        """
        这个模块是为了清洗因为周期性边界条件所带来的原子坐标跃迁
        :param coord_tmp: 在某一个时间步下体系原子的坐标
        :param coord_ini: 体系最初始的原子坐标
        :param cell_tmp: 晶格常数矩阵
        :return: 消除因周期性边界条件带来的跃迁后的后的原子坐标矩阵coord_ini
        """
        cell_tmp_inv = np.linalg.inv(cell_tmp)  # cell_tmp求逆
        coord_diff = coord_tmp - coord_ini  # 判断依据，利用与初始坐标的差值来判断
        coord_diff_frac = np.dot(coord_diff, cell_tmp_inv)
        # 这一步是将直角坐标转换为分数坐标，coord_diff_frac * cell_tmp=coord_diff，cell_tmp_inv是cell_tmp的逆矩阵
        coord_tmp_frac = np.dot(coord_ini, cell_tmp_inv)  # fractional coordinates
        # 这一步是将直角坐标转换为分数坐标，coord_tmp_frac * cell_tmp=coord_tmp，cell_tmp_inv是cell_tmp的逆矩阵
        for k in range(3):
            index1 = np.where(coord_diff_frac[:, k] > 0.5)  # periodic boundary condition
            coord_tmp_frac[index1, k] += 1.0
            index2 = np.where(coord_diff_frac[:, k] < -0.5)
            coord_tmp_frac[index2, k] -= 1.0
        coord_wash = np.dot(coord_tmp_frac, cell_tmp)  # 将分数坐标再转化回来
        return coord_wash

    def read_cell_coord(self, f):
        """
        输出某一时间步的晶格常数和消除周期性影响的原子坐标信息
        :return:晶格常数矩阵cell_tmp和原子位置矩阵coord_wash
        """
        #num_files, type_index, coord_ini = self.para_int()
        skip_step = 5  # 跳过前5行文字描述,以便读取晶格常数信息
        for j in range(skip_step):
            f.readline()
        "下面的代码：读取晶格常数的信息"
        cell_tmp = self.__readcell(f)
        "下面的代码：读取原子坐标的信息"
        f.readline()  # 跳过ITEM: ATOMS id type x y z这一行
        type_index, coord_tmp, id_inex = self.__readatoms(f, self.Number_atoms)
        coord_wash = self.__data_wash(self.coord_ini, coord_tmp, cell_tmp)
        return cell_tmp, coord_wash


    def __core_inside(self, criteria, final_coord, initial_coord, index, num, cell):
        """
        MXS也不知道这个函数究竟是什么了
        :param criteria: 最近邻判断标准
        :param final_coord: 原子坐标矩阵输入
        :param initial_coord: 原子坐标矩阵输入
        :param index: 原子的id矩阵
        :param num: 内部参数
        :param cell: 晶格常数矩阵
        :return: 通过最近邻标准筛选后的
        """
        inside_final = np.full((initial_coord.shape[0], self.dim), final_coord[num, :])
        "#建立一个所有行完全一致，列按照final_coord[num, :]排布的矩阵inside_final"
        inside_wash = self.__data_wash(inside_final, initial_coord, cell)
        "#通过inside_final - initial_coord作为判据调用函数data_wash来清洗应周期性条件而产生的原子坐标跃迁，data_wash_switch=1, 输出的是清洗过的initial_coord"
        inside = inside_final - inside_wash
        "#inside代表两个不同类型原子坐标的差值，也就是距离"
        bond_where = np.where(abs(inside) < criteria, 1, 0)
        "#bond_where代表寻找矩阵inside中每个绝对值小于criteria的元素，真则赋值为1，假则赋值为0"
        bond_mask = (np.sum(bond_where, axis=1) == self.dim)
        "#产生模板矩阵bond_mask，当按行排列三个元素都为1时，触发==self.dim，赋值true，此时模板矩阵就会保留相应元素"

        return index[bond_mask], np.sum(bond_mask), inside_wash[bond_mask, :], inside[bond_mask, :]


    def insideMain(self, N, cell, final_coord, initial_coord, index, charge, limit_cord):
        """
        判断离子位移和极化的主要函数，也是这个代码的核心所在
        :param N: 原子的总数目
        :param cell: 晶格常数矩阵
        :param final_coord: 原子坐标矩阵输入
        :param initial_coord: 原子坐标矩阵输入
        :param index:原子的id
        :param charge: 有效电荷，计算极化专用。如果计算离子位移，则有效电荷必须设置为1
        :param limit_cord: 限定的配位数
        :return: 1."displament"输出离子位移矩阵，"Polar"输出某位原子的坐标和有效电荷的乘积；2.自适应配位数（由程序自动判断的配位数）3.（1.）的平均值
        """
        record = np.zeros((N, self.dim))  # 用于记录离子位移
        record_id = np.zeros((N, limit_cord))
        record_conumber = np.zeros(N)  # 用于记录自产生配位数
        def takeZero(elem):
            """这个函数是为了选取矩阵的第0列"""
            return elem[0]
        test_array = [[] for i in range(N)]
        print(len(test_array))
        for num in range(N):
            id, record_conumber[num], inside_wash, inside = self.__core_inside(self.criteria, final_coord, initial_coord, index, num, cell)
            "#np.sum(bond_where, axis=1)代表对bond_where按列求和，如果和值为3就代表xyz三个方向都小于criteria，对应的点就是最近邻点，并得到N*1维的模板矩阵"
            #np.save(f"inside_wash_{num}.npy", inside_wash)
            test_array[num] = inside_wash.copy()
            fin_bond = np.sum(charge[num] * inside_wash, axis=0) / record_conumber[num]
            "#inside_wash[bond_mask, :]代表只保留利用模板矩阵清洗过的initial_coord的数据，此时数据的个数就是配位数，将这些数据按列相加，并乘有效电荷"

            if self.switch == 1:
                if record_conumber[num] < limit_cord:
                    criteria_ = self.criteria
                    while record_conumber[num] < limit_cord:
                        criteria_ = criteria_+0.15
                        id, record_conumber[num], inside_wash, inside = self.__core_inside(criteria_, final_coord, initial_coord, index, num, cell)
                        fin_bond = np.sum(charge[num] * inside_wash, axis=0) / record_conumber[num]
                        #np.save(f"inside_wash_{num}.npy", inside_wash)
                        test_array[num] = inside_wash.copy()


                if record_conumber[num] > limit_cord:
                    insideX = list(np.column_stack((np.sum(inside * inside, axis=1), id, inside_wash)))
                    insideX.sort(key=takeZero)
                    fin = np.array(insideX)[:limit_cord, :]
                    id = fin[:, 1]
                    fin_inside_wash = fin[:, 2:]
                    #np.save(f"inside_wash_{num}.npy", fin_inside_wash)
                    test_array[num] = fin_inside_wash.copy()
                    fin_bond = np.sum(charge[num] * fin_inside_wash, axis=0) / limit_cord
                    if num == N-1:
                        print(fin_bond)
                    record_conumber[num] = limit_cord
            try:
                record_id[num, :] = id
                
            except: self.error = "In fact, the coordination number is incorrect"
            record[num, :] = fin_bond
        return record, record_conumber, record_id


class AVG_disp_polar(__Core__):
    def ALLcoord_avg(self):
        "计算timestep_selected时间步内的平均晶格常数和平均坐标位置"
        cell_avg = np.zeros((3, 3))
        coord_avg = np.zeros((self.Number_atoms, 3))

        with open(self.input_path, 'r') as f:
            for i in range(self.num_files):  # 不同时间步总数目
                cell_tmp, coord_wash = self.read_cell_coord(f)
                if i >= (self.num_files - self.timestep_selected):
                    cell_avg += cell_tmp
                    coord_avg += coord_wash
            cell_avg = cell_avg / self.timestep_selected
            coord_avg = coord_avg / self.timestep_selected
            f.close()
        return cell_avg, coord_avg

    def get_atomsdisplacment(self, cell, coord):
        """
        这个函数用于将原子位置矩阵转化为相对极化位移矩阵
        This program is written by MXS, for sloving the displacemnet of atoms.
        :param cell: 晶格常数矩阵，由read_cell_coord或ALLcoord_avg产生
        :param coord: 原子坐标矩阵，由read_cell_coord或ALLcoord_avg产生
        :return: 原子位移矩阵和自产生配位数
        """
        Acoord = coord[:self.N_A, :]
        Bcoord = coord[self.N_A:(self.N_A + self.N_B), :]
        Xcoord = coord[(self.N_A + self.N_B):, :]
        Xindex = self.id[(self.N_A + self.N_B):]
        # 原子位置坐标
        record_AX, record_conumber_AX, id_AX = self.insideMain(self.N_A, cell, Acoord, Xcoord, Xindex, np.ones(self.N_A), 12)
        record_BX, record_conumber_BX, id_BX = self.insideMain(self.N_B, cell, Bcoord, Xcoord, Xindex, np.ones(self.N_B), 6)
        "#注意，计算极化时charge必须为1"
        disp_A, disp_B = Acoord - record_AX, Bcoord - record_BX
        print(Acoord[-1],record_AX[-1],disp_A[-1])
        disp_fin = np.vstack((disp_A, disp_B))
        corrd_number_fin = np.append(record_conumber_AX, record_conumber_BX)
        record_avg_fin = np.vstack((np.sum(disp_A, axis=0)/self.N_A, np.sum(disp_B, axis=0)/self.N_B))
        return disp_fin, corrd_number_fin, record_avg_fin, id_AX, id_BX, self.error

    def get_totPolar_cell(self, cell, coord):
        """
        get_Polar_cell用于将原子坐标转化为单位晶胞内的极化强度（C/m^2）
        :param cell:
        :param coord:
        :return:
        """
        Acoord = coord[:self.N_A, :]
        Aindex = self.id[:self.N_A]
        Bcoord = coord[self.N_A:(self.N_A + self.N_B), :]
        Xcoord = coord[(self.N_A + self.N_B):, :]
        Xindex = self.id[(self.N_A + self.N_B):]
        "#原子位置坐标"
        charge_A = [self.charge_A[clf] for clf in self.type_index[:self.N_A]]
        charge_B = [self.charge_B[clf] for clf in self.type_index[self.N_A:(self.N_A + self.N_B)]]
        charge_X = np.full(self.N_A, self.charge_X)
        "#离子有效电荷"
        Volum = np.cross(cell[0, :], cell[1, :]) * cell[2, :]
        "#晶胞体积a×b·c，np.cross()是两个向量的叉乘"
        PolarB_A, record_conumberB_A, record_id_BA = self.insideMain(self.N_B, cell, Bcoord, Acoord, Aindex, charge_A, 8)
        PolarB_X, record_conumberB_X, record_id_BX = self.insideMain(self.N_B, cell, Bcoord, Xcoord, Xindex, charge_X, 6)
        PolarB = [Q_B*C_B for Q_B, C_B in zip(charge_B, Bcoord)]

        Polar_cell = (PolarB_A + 3 * PolarB_X + PolarB) * self.cellAll * self.ff / Volum[2]
        return Polar_cell

    def Polar_eachAtoms(self, cell, dispcoord):
        """
        Polar_eachAtoms函数用于计算每个A位和B位分别产生的极化，计算方法是离子位移乘有效电荷
        :param cell:晶胞体积
        :param dispcoord:离子位移
        :return:A位和B位分别产生的极化和平均极化
        """
        Adispcoord = dispcoord[:self.N_A, :]
        Bdispcoord = dispcoord[self.N_A:(self.N_A + self.N_B), :]

        charge_A = [self.charge_A[clf] for clf in self.type_index[:self.N_A]]
        charge_B = [self.charge_B[clf] for clf in self.type_index[self.N_A:(self.N_A + self.N_B)]]

        Volum = np.cross(cell[0, :], cell[1, :]) * cell[2, :]
        "#晶胞体积a×b·c，np.cross()是两个向量的叉乘"

        PolarA = np.array([Q_A * D_A for Q_A, D_A in zip(charge_A, Adispcoord)]) * self.cellAll * self.ff / Volum[2]
        PolarB = np.array([Q_B * D_B for Q_B, D_B in zip(charge_B, Bdispcoord)]) * self.cellAll * self.ff / Volum[2]

        avg_PolarA, avg_PolarB = np.sum(PolarA, axis=0)/self.N_A, np.sum(PolarB, axis=0)/self.N_B
        return PolarA, PolarB, avg_PolarA, avg_PolarB, avg_PolarA + avg_PolarB


class Output:
    def __init__(self, Number_atoms, Number_cell):
        self.Number_atoms = Number_atoms
        self.cellAll = Number_cell[0] * Number_cell[1] * Number_cell[2]

    def output_avg_cell_coord(self, output_path, cell_avg, coord_avg):
        outfile = open(output_path, 'w')
        print("%22s %22s %22s" % ("cell_x", "cell_y", "cell_z"), file=outfile)
        for i in range(3):
            print("%25.16f %22.16f %22.16f" % (cell_avg[i, 0], cell_avg[i, 1], cell_avg[i, 2]), file=outfile)
        print("%22s %22s %22s" % ("coord_x", "coord_y", "coord_z"), file=outfile)
        for i in range(self.Number_atoms):
            print("%22.16f %22.16f %22.16f" % (coord_avg[i, 0], coord_avg[i, 1], coord_avg[i, 2]), file=outfile)
        outfile.close()

    def output_disp(self, output_path, disp_fin, corrd_number_fin, record_avg_fin, error):
        outfile = open(output_path, 'w')
        print(error, file=outfile)
        print("%22s %22s %22s" % ("avgdisp_x", "avgdisp_y", "avgdisp_z"), file=outfile)
        for j in range(2):
            print("%22.16f %22.16f %22.16f" % (record_avg_fin[j, 0], record_avg_fin[j, 1], record_avg_fin[j, 2]), file=outfile)
        print("%22s %22s %22s %22s" % ("displacement_x", "displacement_y", "displacement_z", "coordination number"), file=outfile)
        for i in range(len(corrd_number_fin)):
            print("%22.16f %22.16f %22.16f %22.1f" % (disp_fin[i, 0], disp_fin[i, 1], disp_fin[i, 2], corrd_number_fin[i]), file=outfile)
        outfile.close()

    def output_id(self, output_path, id_AX, id_BX):
        """用于输出配位数和原子id"""
        outfile = open(output_path, 'w')
        for i in range(id_AX.shape[0]):
            print(id_AX[i, :].astype(int), file=outfile)
        for i in range(id_BX.shape[0]):
            print(id_BX[i, :].astype(int), file=outfile)
        outfile.close()

    def output_polar_cell(self, output_path, Polar_cell):
        outfile = open(output_path, 'w')
        print("%22s %22s %22s" % ("Polar_x", "Polar_y", "Polar_z"), file=outfile)
        for i in range(self.cellAll):
            print("%22.16f %22.16f %22.16f" % (Polar_cell[i, 0], Polar_cell[i, 1], Polar_cell[i, 2]), file=outfile)
        outfile.close()

    def output_Polar_eachAtoms(self, output_path, PolarA, PolarB, avg_PolarA, avg_PolarB, avg_AllPolar):
        outfile = open(output_path, 'w')
        print("%22s %22s %22s" % ("avgPolarA_x", "avgPolarA_y", "avgPolarA_z"), file=outfile)
        print("%22.16f %22.16f %22.16f" % (avg_PolarA[0], avg_PolarA[1], avg_PolarA[2]), file=outfile)
        print("%22s %22s %22s" % ("avgPolarB_x", "avgPolarB_y", "avgPolarB_z"), file=outfile)
        print("%22.16f %22.16f %22.16f" % (avg_PolarB[0], avg_PolarB[1], avg_PolarB[2]), file=outfile)
        print("%22s %22s %22s" % ("avgPolarAll_x", "avgPolarAll_y", "avgPolarAll_z"), file=outfile)
        print("%22.16f %22.16f %22.16f" % (avg_AllPolar[0], avg_AllPolar[1], avg_AllPolar[2]), file=outfile)
        print("%22s %22s %22s" % ("Polar_x", "Polar_y", "Polar_z"), file=outfile)
        print("The polarization of A-site ion and B-site ion respectively", file=outfile)
        print("%22s %22s %22s" % ("Polar_x", "Polar_y", "Polar_z"), file=outfile)
        for i in range(self.cellAll):
            print("%22.16f %22.16f %22.16f" % (PolarA[i, 0], PolarA[i, 1], PolarA[i, 2]), file=outfile)
        for i in range(self.cellAll):
            print("%22.16f %22.16f %22.16f" % (PolarB[i, 0], PolarB[i, 1], PolarB[i, 2]), file=outfile)
        outfile.close()

