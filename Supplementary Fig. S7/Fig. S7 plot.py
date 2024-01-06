import matplotlib.pyplot as plt
x = ['Sample 151507', 'Sample 151508', 'Sample 151509', 'Sample 151510', 'Sample 151669', 'Sample 151670',
     'Sample 151671', 'Sample 151672', 'Sample 151673', 'Sample 151674', 'Sample 151675', 'Sample 151676']
y = [0.5966, 0.5451, 0.5638, 0.5296, 0.7405, 0.7807, 0.8243, 0.7812, 0.6228, 0.5695, 0.5486, 0.5169]
z = [0.5054, 0.4763, 0.4391, 0.4395, 0.6223, 0.4733, 0.5698, 0.7670, 0.5275, 0.5176, 0.5450, 0.4755]



# 增加一个固定维度，长度与上述数据一样
fix_value = []
# 求出数据y和z的最大值，取其1/4的值作为固定长度
value_max = max(max(y), max(z))
fix_temp = value_max / 4
# fix_temp = value_max
for i in range(len(x)):
 fix_value.append(fix_temp)
print(fix_value)
# 将y，z其中一个维度变为负值，我们选择z
z_ne = [-i for i in z]
print(z_ne)

'''
# 设置中文显示为微软雅黑
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
'''
# 设置图表默认字体为12
plt.rcParams['font.size'] = 11
# 设置画布大小
plt.figure(figsize=(14, 8))
# 画条形图,设置颜色和柱宽，将fix_value,y,z_ne依次画进图中
plt.barh(x, fix_value, color='w', height=0.5)
plt.barh(x, y, left=fix_value, color='#82B0D2', label='Augment', height=0.5)
plt.barh(x, z_ne, color='#BEB8DC', height=0.5, label='No augment')

# 添加数据标签，将fix_value的x标签显示出来，y和z_ne的数据标签显示出来
for a, b in zip(x, fix_value):
 plt.text(b / 2, a, str(a), ha='center', va='center', fontsize=12)
 # plt.text(b/2, a, '%s' % str(a), ha='center', va='center', fontsize=12)
for a, b in zip(x, y):
 plt.text(b + fix_temp + value_max / 20, a, str(b), ha='center', va='center')
 # plt.text(b + fix_temp + value_max / 20, a, '%d' % int(b), ha='center', va='center')
for a, b in zip(x, z):
 plt.text(-b - value_max / 20, a, str(b), ha='center', va='center')
 # plt.text(-b - value_max / 20, a, '%d' % int(b), ha='center', va='center')

# 坐标轴刻度不显示
plt.xticks([])
plt.yticks([])
# 添加图例，自定义位置
# plt.legend(bbox_to_anchor=(-0.02, 0.5), frameon=False, fontsize=12)
plt.legend(frameon=False, fontsize=12)
# 添加标题，并设置字体
# plt.title(label='比比谁家玩具多呀', fontsize=12, fontweight='bold')
# 设置绘图区域边框不可见
ax = plt.gca()
ax.set_axisbelow(True)
# 设置绘图区域边框不可见
[ax.spines[loc_axis].set_visible(False) for loc_axis in ['bottom', 'top', 'right', 'left']]
# 使布局更合理
plt.tight_layout()
plt.show()
plt.savefig('./augment test/bar plot.pdf', dpi=300)






