Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(24, 29, 36)
$muted = [System.Drawing.Color]::FromArgb(90, 98, 108)
$grid = [System.Drawing.Color]::FromArgb(210, 215, 222)
$point = [System.Drawing.Color]::FromArgb(115, 98, 98, 98)
$line = [System.Drawing.Color]::FromArgb(20, 24, 30)

function New-Font($size, $style = [System.Drawing.FontStyle]::Regular) {
  return [System.Drawing.Font]::new($fontFamily, [single]$size, $style)
}

function New-Brush($color) {
  return [System.Drawing.SolidBrush]::new($color)
}

function New-Pen($color, $width = 1) {
  return [System.Drawing.Pen]::new($color, [single]$width)
}

function Draw-Text($g, $text, $x, $y, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $g.DrawString($text, $font, $brush, [single]$x, [single]$y)
  $font.Dispose()
  $brush.Dispose()
}

function Draw-Centered($g, $text, $x, $y, $w, $h, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $fmt = [System.Drawing.StringFormat]::new()
  $fmt.Alignment = [System.Drawing.StringAlignment]::Center
  $fmt.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $font, $brush, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $fmt)
  $fmt.Dispose()
  $font.Dispose()
  $brush.Dispose()
}

function Nice-Range($values) {
  $min = ($values | Measure-Object -Minimum).Minimum
  $max = ($values | Measure-Object -Maximum).Maximum
  if ($min -eq $max) {
    return @([double]$min - 1, [double]$max + 1)
  }
  $pad = ([double]$max - [double]$min) * 0.08
  return @(([double]$min - $pad), ([double]$max + $pad))
}

function Scale-X($x, $xmin, $xmax, $left, $right) {
  return $left + (([double]$x - [double]$xmin) / ([double]$xmax - [double]$xmin)) * ($right - $left)
}

function Scale-Y($y, $ymin, $ymax, $top, $bottom) {
  return $bottom - (([double]$y - [double]$ymin) / ([double]$ymax - [double]$ymin)) * ($bottom - $top)
}

function Normal-Cdf($z) {
  $sign = 1.0
  if ($z -lt 0) { $sign = -1.0 }
  $x = [Math]::Abs([double]$z) / [Math]::Sqrt(2.0)
  $t = 1.0 / (1.0 + 0.3275911 * $x)
  $a1 = 0.254829592
  $a2 = -0.284496736
  $a3 = 1.421413741
  $a4 = -1.453152027
  $a5 = 1.061405429
  $erf = 1.0 - (((((($a5 * $t + $a4) * $t) + $a3) * $t + $a2) * $t + $a1) * $t * [Math]::Exp(-$x * $x))
  return 0.5 * (1.0 + $sign * $erf)
}

function Format-P($p) {
  $v = [double]$p
  if ($v -lt 0.001) { return "p < 0.001" }
  return ("p = {0:0.003}" -f $v)
}

function Get-PanelStats($rows, $metric) {
  $n = @($rows).Count
  $xs = @($rows | ForEach-Object { [double]$_.mean_reflectance })
  $ys = @($rows | ForEach-Object { [double]$_.$metric })
  $xMean = ($xs | Measure-Object -Average).Average
  $yMean = ($ys | Measure-Object -Average).Average
  $sxx = 0.0
  $syy = 0.0
  $sxy = 0.0
  for ($i = 0; $i -lt $n; $i++) {
    $dx = $xs[$i] - $xMean
    $dy = $ys[$i] - $yMean
    $sxx += $dx * $dx
    $syy += $dy * $dy
    $sxy += $dx * $dy
  }
  $slope = $sxy / $sxx
  $intercept = $yMean - $slope * $xMean
  $r = $sxy / [Math]::Sqrt($sxx * $syy)
  $r2 = $r * $r
  $t = [Math]::Abs($r) * [Math]::Sqrt(($n - 2) / (1 - $r2))
  $pApprox = 2 * (1 - (Normal-Cdf $t))
  $pApprox = [Math]::Max(0.0, [Math]::Min(1.0, $pApprox))
  return [pscustomobject]@{
    n = $n
    pearson_r = $r
    r_squared = $r2
    p_value = $pApprox
    slope = $slope
    intercept = $intercept
  }
}

function Draw-Panel($g, $rows, $metric, $panelTitle, $x, $y, $w, $h, $showYLabel, $showXLabel) {
  $left = $x + 96
  $right = $x + $w - 36
  $top = $y + 72
  $bottom = $y + $h - 84

  $xs = @($rows | ForEach-Object { [double]$_.mean_reflectance })
  $ys = @($rows | ForEach-Object { [double]$_.$metric })
  $xr = Nice-Range $xs
  $yr = Nice-Range $ys
  $xmin = $xr[0]; $xmax = $xr[1]
  $ymin = $yr[0]; $ymax = $yr[1]

  Draw-Centered $g $panelTitle $x ($y + 12) $w 42 23 $ink ([System.Drawing.FontStyle]::Bold)

  $gridPen = New-Pen $grid 1
  $axisPen = New-Pen $ink 2
  for ($i = 0; $i -le 4; $i++) {
    $tx = $left + (($right - $left) * $i / 4)
    $ty = $top + (($bottom - $top) * $i / 4)
    $g.DrawLine($gridPen, $tx, $top, $tx, $bottom)
    $g.DrawLine($gridPen, $left, $ty, $right, $ty)
  }
  $g.DrawLine($axisPen, $left, $bottom, $right, $bottom)
  $g.DrawLine($axisPen, $left, $top, $left, $bottom)

  for ($i = 0; $i -le 4; $i++) {
    $xv = $xmin + (($xmax - $xmin) * $i / 4)
    $tx = $left + (($right - $left) * $i / 4)
    Draw-Centered $g ("{0:0.000}" -f $xv) ($tx - 48) ($bottom + 12) 96 28 15 $muted
    $yv = $ymin + (($ymax - $ymin) * $i / 4)
    $ty = $bottom - (($bottom - $top) * $i / 4)
    Draw-Text $g ("{0:0.0}" -f $yv) ($left - 88) ($ty - 12) 15 $muted
  }

  $ptBrush = New-Brush $point
  foreach ($row in $rows) {
    $px = Scale-X $row.mean_reflectance $xmin $xmax $left $right
    $py = Scale-Y $row.$metric $ymin $ymax $top $bottom
    $g.FillEllipse($ptBrush, $px - 4, $py - 4, 8, 8)
  }
  $ptBrush.Dispose()

  $stat = Get-PanelStats $rows $metric
  $slope = [double]$stat.slope
  $intercept = [double]$stat.intercept
  $y1 = $intercept + $slope * $xmin
  $y2 = $intercept + $slope * $xmax
  $x1p = Scale-X $xmin $xmin $xmax $left $right
  $x2p = Scale-X $xmax $xmin $xmax $left $right
  $y1p = Scale-Y $y1 $ymin $ymax $top $bottom
  $y2p = Scale-Y $y2 $ymin $ymax $top $bottom
  $linePen = New-Pen $line 3
  $g.DrawLine($linePen, $x1p, $y1p, $x2p, $y2p)
  $linePen.Dispose()
  $statText = ("R2 = {0:0.000}; {1}" -f [double]$stat.r_squared, (Format-P $stat.p_value))
  Draw-Text $g $statText ($left + 14) ($top + 14) 18 $ink ([System.Drawing.FontStyle]::Bold)

  if ($showXLabel) {
    Draw-Centered $g "Overall brightness (mean reflectance)" ($left + 40) ($bottom + 42) (($right - $left) - 80) 36 18 $ink
  }
  if ($showYLabel) {
    $state = $g.Save()
    $g.TranslateTransform($x + 18, $top + (($bottom - $top) / 2) + 98)
    $g.RotateTransform(-90)
    Draw-Centered $g $panelTitle 0 0 196 34 17 $ink
    $g.Restore($state)
  }

  $gridPen.Dispose()
  $axisPen.Dispose()
}

$dataset = Import-Csv (Join-Path $root "reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_dataset.csv")
$reflectance = Import-Csv (Join-Path $root "reports/tables/pc1_mean_reflectance_correlation/50m_regular_PCA_pc_mean_reflectance_quadrat_means.csv")

$reflectanceByKey = @{}
foreach ($r in $reflectance) {
  $reflectanceByKey["$($r.scale)|$($r.quad_id)"] = $r
}

$rows = foreach ($d in $dataset) {
  if ($d.scale -ne "50m") { continue }
  $key = "$($d.scale)|$($d.quad_id)"
  if (-not $reflectanceByKey.ContainsKey($key)) { continue }
  $r = $reflectanceByKey[$key]
  $needed = @(
    $r.mean_reflectance,
    $d.spec_pca_mean,
    $d.spec_spca_mean,
    $d.spec_rao_q,
    $d.spec_spca_rao,
    $d.spec_alpha,
    $d.spec_spca_alpha
  )
  if ($needed | Where-Object { [string]::IsNullOrWhiteSpace([string]$_) -or [string]$_ -eq "NA" }) { continue }
  [pscustomobject]@{
    quad_id = $d.quad_id
    scale = $d.scale
    mean_reflectance = [double]$r.mean_reflectance
    spec_pca_mean = [double]$d.spec_pca_mean
    spec_spca_mean = [double]$d.spec_spca_mean
    spec_rao_q = [double]$d.spec_rao_q
    spec_spca_rao = [double]$d.spec_spca_rao
    spec_alpha = [double]$d.spec_alpha
    spec_spca_alpha = [double]$d.spec_spca_alpha
  }
}

$width = 2200
$height = 2500
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "PCA spectral heterogeneity metrics and overall brightness" 110 64 36 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "50 m quadrats; overall brightness is mean reflectance across retained wavelengths" 110 112 22 $muted
Draw-Centered $g "Raw PCA" 230 170 760 44 28 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Centered $g "Standardized PCA" 1190 170 760 44 28 $ink ([System.Drawing.FontStyle]::Bold)

$panelW = 900
$panelH = 680
$xRaw = 120
$xStd = 1080
$y0 = 230
$gapY = 70

$panels = @(
  @{ RawMetric = "spec_pca_mean"; StdMetric = "spec_spca_mean"; Title = "mean distance" },
  @{ RawMetric = "spec_rao_q"; StdMetric = "spec_spca_rao"; Title = "Rao's Q" },
  @{ RawMetric = "spec_alpha"; StdMetric = "spec_spca_alpha"; Title = "alpha hull" }
)

for ($i = 0; $i -lt $panels.Count; $i++) {
  $y = $y0 + $i * ($panelH + $gapY)
  $showX = ($i -eq $panels.Count - 1)
  Draw-Panel $g $rows $panels[$i].RawMetric $panels[$i].Title $xRaw $y $panelW $panelH $true $showX
  Draw-Panel $g $rows $panels[$i].StdMetric $panels[$i].Title $xStd $y $panelW $panelH $false $showX
}

$outPath = Join-Path $outDir "09_raw_vs_standardized_pca_metric_brightness_scatter_50m.png"
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Host $outPath
