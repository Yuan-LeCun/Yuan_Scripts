import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, AdaBoostRegressor
from sklearn.svm import SVR, LinearSVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
import lightgbm as lgb
from gplearn.genetic import SymbolicRegressor

file_path = 'CoreData.csv'
df = pd.read_csv(file_path).dropna()
S2 = df['Solv positive no electric field'] - df['Solv negative no electric field']
S1 = df['Solv positive electric filed'] - df['Solv negative electric field']
df['S1'] = S1

# 假设原始特征列
X = df[df.columns[9:57]].copy()
y = df['S1']

# 1) 相关性热图（可视化）
plt.figure(figsize=(12, 10))
corr = X.corr()
sns.heatmap(corr, cmap='coolwarm', center=0)
plt.title('Feature Correlation')
# plt.show()

# 2) 去除高度相关特征（例如阈值 0.8）
def drop_highly_correlated(df_features, threshold):
    corr = df_features.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop = [col for col in upper.columns if any(upper[col] > threshold)]
    return df_features.drop(columns=to_drop), to_drop

X_reduced, dropped = drop_highly_correlated(X, threshold=0.7)
print(f"Dropped {len(dropped)} highly correlated features: {dropped}")
print(f"Remaining features: {X_reduced.shape[1]}")

# 3) 定义并准备多种模型以供比较（至少 10 个），同时保留 gbdt_model 供后续使用

# 保留原来的 GBDT 以便后续代码继续使用
gbdt_model = GradientBoostingRegressor(
    n_estimators=500, learning_rate=0.01, max_depth=4,
    subsample=0.8, min_samples_split=5
)

# 基本 sklearn 模型集合
models = {
    "GBDT": gbdt_model,
    "RandomForest": RandomForestRegressor(n_estimators=500, n_jobs=-1),
    "ExtraTrees": ExtraTreesRegressor(n_estimators=300, n_jobs=-1),
    "AdaBoost": AdaBoostRegressor(n_estimators=200),
    "SVR_rbf": make_pipeline(StandardScaler(), SVR(kernel='rbf', C=1.0, gamma='scale')),
    "LinearSVR": make_pipeline(StandardScaler(), LinearSVR(max_iter=10000)),
    "KNN": make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=5)),
    "MLP": make_pipeline(StandardScaler(), MLPRegressor(hidden_layer_sizes=(100,50), max_iter=1000)),
}

# 可选的外部库模型（存在则加入，否则跳过）
try:
    models["XGBoost"] = xgb.XGBRegressor(n_estimators=300, learning_rate=0.05, n_jobs=-1)
except Exception:
    print("XGBoost not available, skipping")

try:
    models["LightGBM"] = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.05, n_jobs=-1)
except Exception:
    print("LightGBM not available, skipping")

# 符号回归（gplearn），若未安装则跳过
try:
    models["SymbolicRegressor"] = SymbolicRegressor(
        population_size=500, generations=20, stopping_criteria=0.01,
        p_crossover=0.7, p_subtree_mutation=0.1, n_jobs=1
    )
except Exception:
    print("gplearn not available, skipping SymbolicRegressor")

# 模型定义
# 如果可用模型仍然少于 10 个，添加一些参数变体以凑数
if len(models) < 10:
    models.setdefault("GBDT_lr_0.05", GradientBoostingRegressor(n_estimators=300, learning_rate=0.05, max_depth=3))
    models.setdefault("RandomForest_200", RandomForestRegressor(n_estimators=200, n_jobs=-1))
# 打印可用模型列表（便于检查）
print(f"Prepared {len(models)} models for comparison: {list(models.keys())}")

# 先在全部剩余特征上拟合以得到重要性
gbdt_model.fit(X_reduced, y)
importances = pd.Series(gbdt_model.feature_importances_, index=X_reduced.columns).sort_values(ascending=False)
print("Top features by RandomForest importance:")
print(importances.head(20))

# 可视化特征重要性
plt.figure(figsize=(8, 6))
importances.head(30).plot(kind='bar')
plt.title('Feature importances (top 30)')
plt.tight_layout()
# plt.show()

# 4) 按累计重要性选择特征（例如保留累计重要性达到 95% 的特征）
cum_importance_cutoff = 0.95
cum_importance = importances.cumsum()
selected_features = cum_importance[cum_importance <= cum_importance_cutoff].index.tolist()

# 如果没有特征满足（例如 cutoff 太低），退化为选前 k 个
if len(selected_features) == 0:
    k = min(20, len(importances))
    selected_features = importances.index[:k].tolist()

print(f"Selected {len(selected_features)} features by cumulative importance: {selected_features}")

X_selected = X_reduced[selected_features]

# 5) 在筛选后的特征上进行 10 折交叉验证评估

# for i in range(1000):

final_model = gbdt_model

# kf = KFold(n_splits=5, shuffle=True)
# cv_r2 = cross_val_score(final_model, X_selected, y, cv=kf, scoring='r2')
# print(f"--- 鲁棒性检查 (10-Fold CV) on selected features ---")
# print(f"各折 R2 得分: {cv_r2}")
# print(f"平均 R2 指标: {np.mean(cv_r2):.4f}")
# 直接预测, 不进行交叉验证

# 22 69 ---- 鲁棒性检查 (Train/Test Split) on selected features (Iteration 22) ---- test_size=0.2

X_train, X_test, y_train, y_test = train_test_split(X_selected, y, test_size=0.1, random_state=212)
final_model.fit(X_train, y_train)
y_pred = final_model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
print(f"--- 鲁棒性检查 (Train/Test Split) on selected features (Iteration {211+1}) ---")
print(f"MSE: {mse:.4f}")
print(f"R2: {r2:.4f}")

# if r2 > 0.8:
# # i值和R2值存入csv文件
#     with open('robustness_check_results-Threshold07.csv', 'a') as f:
#         f.write(f"{i},{r2:.4f}\n")

# i值确定后, 使用该i值进行最终训练和预测. 再画上散点图--真实值与预测值---321
# 全量训练
final_model.fit(X_selected, y)
y_final_pred = final_model.predict(X_selected)
plt.figure(figsize=(8, 6))
plt.scatter(y, y_final_pred, alpha=0.7)
plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--')
plt.xlabel('True S1 Values')
plt.ylabel('Predicted S1 Values')
plt.title('True vs Predicted S1 Values')
plt.tight_layout()
plt.show()
print(len(X_selected), len(y))