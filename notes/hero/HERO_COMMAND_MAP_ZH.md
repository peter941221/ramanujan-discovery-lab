# Hero 命令地图

状态日期: `2026-03-20`

## 目的

这份短文档回答一个很实用的问题:

```text
当我们围着 `cb60fd71d1d7` 工作时
每个 CLI 命令到底负责什么
它更像哪一种研究动作
什么时候该用它
什么时候不该用它
```

## 一张总图

```text
Hero 工作流
├─ `discover`
│  └─ 海选候选
├─ `verify`
│  └─ 做高精度复核
├─ `analyze`
│  └─ 做局部对比体检
├─ `research`
│  └─ 做重型专项化验
├─ `identify`
│  └─ 做身份鉴定 / 溯源
├─ `tail-note`
│  └─ 专看尾部家族像谁
├─ `tail-operator-note`
│  └─ 专看尾部是否藏着算子规律
├─ `formalize`
│  └─ 把当前结论翻成证明工作清单
├─ `report`
│  └─ 产出公开摘要
└─ `site`
   └─ 渲染成静态站点
```

## 基础词树

```text
核心术语
├─ candidate
│  ├─ 中文: 候选对象
│  └─ 类比: 像一块刚从矿里挖出来的石头, 先知道它值得看, 还不知道是不是宝石
├─ benchmark
│  ├─ 中文: 已知基准对象
│  └─ 类比: 博物馆里已经编号入库的标准样本
├─ identify
│  ├─ 中文: 身份鉴定 / 源对象识别
│  └─ 类比: 拿着候选去和已知家族、配方、指纹库逐项比对
├─ tail family
│  ├─ 中文: 尾部家族
│  └─ 类比: 看整条 continued fraction 的“后半段骨架”是不是在重复某种生长规律
├─ operator
│  ├─ 中文: 算子 / 递推作用规则
│  └─ 类比: 不是问“它像谁”, 而是问“它遵守什么运动方程”
└─ formalize
   ├─ 中文: 形式化
   └─ 类比: 把研究笔记改写成法院可采纳的逐条证词
```

## 命令职责

### `discover`

- 作用:
  从模板空间里搜索可能有意思的 continued-fraction 候选。
- 输入:
  搜索深度、预算、数值精度。
- 输出:
  候选列表, 通常是 `results/candidates.jsonl`。
- 什么时候用:
  你还在找新目标。
- 什么时候别先用:
  hero case 已经固定, 现在主要是在解释而不是再海选。

### `verify`

- 作用:
  对 `discover` 找到的候选做更高精度数值复核。
- 输出:
  `results/verified.jsonl`
- 类比:
  `discover` 像初筛, `verify` 像复检。

### `analyze`

- 作用:
  给单个候选写一份局部、可读性强的对比说明。
- 重点:
  更像“症状描述”和“对照观察”, 不追求大规模盒子扫描。
- 适用:
  想快速看 hero case 和最近 benchmark 到底差在哪几项低阶系数。

### `research`

- 作用:
  对单个候选跑更重的结构探查。
- 重点:
  会做较多符号级实验, 包括 product fit、page-43、Bauer-Muir、子序列等。
- 类比:
  这是“重型化验室”。
- 适用:
  你要扩张证据面, 不是只想要一个简短摘要。

### `identify`

- 作用:
  试着把候选写成某个已知对象, 或者至少把“它不像谁”系统地记下来。
- 输出:
  一份身份鉴定笔记, 对 hero case 是
  `notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`
- 它真正问的问题是:

```text
这个 candidate
├─ 是不是某个已知 benchmark
├─ 是不是某个已知 source family 加一个小修正
├─ 是不是 eta-quotient / modular-unit 之类的闭式对象
└─ 如果都不是, 它具体在哪些坐标和盒子里失败
```

- 类比:
  不是“看起来眼熟吗”, 而是“拿指纹、DNA、家谱和历史档案逐项核对”。
- 当前 hero 语境里:
  `identify` 是主力“源对象识别”命令。

### `tail-note`

- 作用:
  不从整个 ratio object 的正面去看, 而是从 exact tail family 的后半段结构去看。
- 关注对象:
  `U_t2`, `U_t3`, `U_t4` 以及 gap-normalized residuals。
- 类比:
  整车难认时, 先看发动机后半段的传动结构像不像某个品牌。
- 适用:
  当普通 `identify` 已经告诉你“正面扫描没认出来”, 但尾部递推结构很有规律时。

### `tail-operator-note`

- 作用:
  在 tail family 上检查低复杂度 q-difference / Mahler / operator 规律。
- 它不主要回答:
  “它像谁?”
- 它主要回答:
  “它是否已经在 obey 某个短小的递推算子?”
- 类比:
  不再查身份证, 而是在测它的运动学方程。
- 适用:
  当我们怀疑最终证明会更像 operator-factorization 或 uniqueness argument。

### `formalize`

- 作用:
  把当前研究状态翻译成 Lean-ready 的形式化准备材料和生成模块。
- 输出:
  hero note:
  `notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md`
- 输出:
  generated Lean file:
  `proofs/Proofs/Generated/Cb60fd71d1d7.lean`
- 它不是:
  最终自动证明机器。
- 它是:
  把“哪些对象、哪些局部引理、哪些 waypoint 已经存在”整理成证明施工图。

### `report`

- 作用:
  汇总当前 verified 结果, 生成对外或对内都更容易扫读的 Markdown 报告。

### `site`

- 作用:
  把结果渲染到静态页面, 供 GitHub Pages 或本地浏览使用。

## Hero 主线推荐顺序

```text
当 hero case 已经固定为 `cb60fd71d1d7`
├─ 日常主线
│  ├─ `analyze`
│  ├─ `research`
│  ├─ `identify`
│  ├─ `tail-note`
│  ├─ `tail-operator-note`
│  └─ `formalize`
└─ 旁线
   ├─ `report`
   └─ `site`
```

- 如果问题是“它像谁”: 先想 `identify`
- 如果问题是“尾部骨架像谁”: 先想 `tail-note`
- 如果问题是“它 obey 什么递推”: 先想 `tail-operator-note`
- 如果问题是“现在能 formalize 到哪”: 先想 `formalize`

## 当前环境下的实用提醒

- `identify` 全量 hero 刷新是长任务, 目前大约几十分钟级别。
- `formalize` 全量 hero 刷新也是长任务, 但通常比 full `identify` 更短。
- `tail-note` 是中等偏重任务。
- `tail-operator-note` 相对更轻。
- 想看 plumbing 是否通, 优先用 `--smoke`。

## 推荐阅读顺序

1. `notes/hero/CB60FD71D1D7_PROGRESS_BOARD.md`
2. `notes/hero/CB60FD71D1D7_AWARD_TRACK_GATE_DASHBOARD.md`
3. `notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`
4. `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
5. `notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
6. `notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md`
