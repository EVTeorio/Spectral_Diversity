Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$tableDir = Join-Path $root "reports/tables/pc1_mean_reflectance_correlation"
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(22, 26, 32)
$muted = [System.Drawing.Color]::FromArgb(83, 90, 99)
$grid = [System.Drawing.Color]::FromArgb(218, 223, 230)
$point = [System.Drawing.Color]::FromArgb(105, 70, 76, 84)
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
  $nums = @($values | ForEach-Object { [double]$_ })
  $min = ($nums | Measure-Object -Minimum).Minimum
  $max = ($nums | Measure-Object -Maximum).Maximum
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

function Get-Regression($rows) {
  $n = [double]$rows.Count
  $sumX = 0.0
  $sumY = 0.0
  foreach ($row in $rows) {
    $sumX += [double]$row.mean_reflectance
    $sumY += [double]$row.pc1
  }
  $meanX = $sumX / $n
  $meanY = $sumY / $n
  $ssXX = 0.0
  $ssXY = 0.0
  foreach ($row in $rows) {
    $dx = [double]$row.mean_reflectance - $meanX
    $dy = [double]$row.pc1 - $meanY
    $ssXX += $dx * $dx
    $ssXY += $dx * $dy
  }
  $slope = $ssXY / $ssXX
  $intercept = $meanY - ($slope * $meanX)
  return @{ Slope = $slope; Intercept = $intercept }
}

function Format-P($p) {
  $value = [double]$p
  if ($value -eq 0 -or $value -lt 0.001) { return "p < 0.001" }
  return ("p = {0:0.000}" -f $value)
}

function Format-R2($r2) {
  return ("R2 = {0:0.000}" -f ([double]$r2))
}

function Import-NumericCsv($path) {
  return @(Import-Csv -Path $path | Where-Object {
    $_.mean_reflectance -ne "" -and $_.mean_reflectance -ne "NA" -and
    $_.pc1 -ne "" -and $_.pc1 -ne "NA"
  })
}

function Get-StatRow($rows, $level) {
  return @($rows | Where-Object { $_.analysis_level -eq $level -and $_.pc_axis -eq "PC1" })[0]
}

function Draw-Panel($g, $panel) {
  $left = $panel.Left
  $top = $panel.Top
  $right = $panel.Right
  $bottom = $panel.Bottom
  $rows = $panel.Rows

  $xRange = Nice-Range @($rows | ForEach-Object { $_.mean_reflectance })
  $yRange = Nice-Range @($rows | ForEach-Object { $_.pc1 })
  $xmin = $xRange[0]
  $xmax = $xRange[1]
  $ymin = $yRange[0]
  $ymax = $yRange[1]

  $axisPen = New-Pen $ink 2.2
  $gridPen = New-Pen $grid 1.1
  $pointBrush = New-Brush $point
  $linePen = New-Pen $line 3

  for ($i = 0; $i -le 4; $i++) {
    $gx = $left + (($right - $left) * $i / 4.0)
    $gy = $top + (($bottom - $top) * $i / 4.0)
    $g.DrawLine($gridPen, [single]$gx, [single]$top, [single]$gx, [single]$bottom)
    $g.DrawLine($gridPen, [single]$left, [single]$gy, [single]$right, [single]$gy)
  }

  foreach ($row in $rows) {
    $px = Scale-X ([double]$row.mean_reflectance) $xmin $xmax $left $right
    $py = Scale-Y ([double]$row.pc1) $ymin $ymax $top $bottom
    $size = if ($panel.Level -eq "pixel") { 4.0 } else { 7.5 }
    $g.FillEllipse($pointBrush, [single]($px - ($size / 2.0)), [single]($py - ($size / 2.0)), [single]$size, [single]$size)
  }

  $regression = Get-Regression $rows
  $lineX1 = $xmin
  $lineX2 = $xmax
  $lineY1 = $regression.Intercept + ($regression.Slope * $lineX1)
  $lineY2 = $regression.Intercept + ($regression.Slope * $lineX2)
  $g.DrawLine(
    $linePen,
    [single](Scale-X $lineX1 $xmin $xmax $left $right),
    [single](Scale-Y $lineY1 $ymin $ymax $top $bottom),
    [single](Scale-X $lineX2 $xmin $xmax $left $right),
    [single](Scale-Y $lineY2 $ymin $ymax $top $bottom)
  )

  $g.DrawLine($axisPen, [single]$left, [single]$bottom, [single]$right, [single]$bottom)
  $g.DrawLine($axisPen, [single]$left, [single]$top, [single]$left, [single]$bottom)

  for ($i = 0; $i -le 4; $i++) {
    $xVal = $xmin + (($xmax - $xmin) * $i / 4.0)
    $yVal = $ymin + (($ymax - $ymin) * $i / 4.0)
    $xPos = $left + (($right - $left) * $i / 4.0)
    $yPos = $bottom - (($bottom - $top) * $i / 4.0)
    Draw-Centered $g ("{0:0.00}" -f $xVal) ($xPos - 42) ($bottom + 10) 84 28 17 $muted
    Draw-Centered $g ("{0:0.0}" -f $yVal) ($left - 86) ($yPos - 14) 74 28 17 $muted
  }

  Draw-Text $g ((Format-R2 $panel.Stat.r_squared) + "; " + (Format-P $panel.Stat.p_value)) ($left + 18) ($top + 16) 22 $ink ([System.Drawing.FontStyle]::Bold)

  $axisPen.Dispose()
  $gridPen.Dispose()
  $pointBrush.Dispose()
  $linePen.Dispose()
}

$rawPixel = Import-NumericCsv (Join-Path $tableDir "50m_regular_PCA_pc_mean_reflectance_pixel_sample.csv")
$rawQuadrat = Import-NumericCsv (Join-Path $tableDir "50m_regular_PCA_pc_mean_reflectance_quadrat_means.csv")
$standardPixel = Import-NumericCsv (Join-Path $tableDir "50m_standardized_PCA_pc_mean_reflectance_pixel_sample.csv")
$standardQuadrat = Import-NumericCsv (Join-Path $tableDir "50m_standardized_PCA_pc_mean_reflectance_quadrat_means.csv")
$rawStats = @(Import-Csv -Path (Join-Path $tableDir "50m_regular_PCA_pc_mean_reflectance_correlation_summary.csv"))
$standardStats = @(Import-Csv -Path (Join-Path $tableDir "50m_standardized_PCA_pc_mean_reflectance_correlation_summary.csv"))

$width = 2400
$height = 1800
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "PC1 and overall brightness" 150 80 44 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "50 m quadrats; overall brightness is mean reflectance across retained wavelengths" 150 140 25 $muted
Draw-Centered $g "Raw PCA" 330 210 760 48 30 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Centered $g "Vector-normalized PCA" 1300 210 760 48 30 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Pixel level" 310 260 27 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Quadrat-mean level" 310 965 27 $ink ([System.Drawing.FontStyle]::Bold)

$panels = @(
  @{ Left = 310; Top = 300; Right = 1090; Bottom = 810; Rows = $rawPixel; Level = "pixel"; Stat = (Get-StatRow $rawStats "pixel_level") },
  @{ Left = 1280; Top = 300; Right = 2060; Bottom = 810; Rows = $standardPixel; Level = "pixel"; Stat = (Get-StatRow $standardStats "pixel_level") },
  @{ Left = 310; Top = 1005; Right = 1090; Bottom = 1515; Rows = $rawQuadrat; Level = "quadrat"; Stat = (Get-StatRow $rawStats "quadrat_mean_level") },
  @{ Left = 1280; Top = 1005; Right = 2060; Bottom = 1515; Rows = $standardQuadrat; Level = "quadrat"; Stat = (Get-StatRow $standardStats "quadrat_mean_level") }
)

foreach ($panel in $panels) {
  Draw-Panel $g $panel
}

Draw-Centered $g "Overall brightness (mean reflectance)" 555 1630 1280 46 27 $ink ([System.Drawing.FontStyle]::Bold)

$font = New-Font 27 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $ink
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("PC1 score", $font, $brush, [System.Drawing.RectangleF]::new(40, 580, 52, 660), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

$outPath = Join-Path $outDir "09_pc1_mean_reflectance_raw_vs_vector_normalized_50m.png"
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Output $outPath
