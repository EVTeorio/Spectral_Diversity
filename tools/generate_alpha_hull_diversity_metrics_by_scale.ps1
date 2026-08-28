Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$tableDir = Join-Path $root "reports/tables/spectral_biodiversity_all_metrics"
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

function Format-P($p) {
  $value = [double]$p
  if ($value -eq 0 -or $value -lt 0.001) { return "p < 0.001" }
  return ("p = {0:0.003}" -f $value)
}

function Format-R($r) {
  return ("r = {0:0.000}" -f ([double]$r))
}

function Format-Tick($value, $metric) {
  if ($metric -in @("sp_shannon", "sp_simpson")) {
    return ("{0:0.00}" -f ([double]$value))
  }
  return ("{0:0}" -f ([double]$value))
}

function Get-StatRow($rows, $scale, $diversityMetric) {
  return @($rows | Where-Object {
    $_.scenario -eq "identity" -and
    $_.scale -eq $scale -and
    $_.spectral_metric -eq "spec_spca_alpha" -and
    $_.biodiversity_metric -eq $diversityMetric
  })[0]
}

function Draw-Panel($g, $panel) {
  $left = $panel.Left
  $top = $panel.Top
  $right = $panel.Right
  $bottom = $panel.Bottom
  $rows = $panel.Rows
  $xMetric = $panel.DiversityMetric

  $xRange = Nice-Range @($rows | ForEach-Object { $_.$xMetric })
  $yRange = Nice-Range @($rows | ForEach-Object { $_.spec_spca_alpha })
  $xmin = $xRange[0]
  $xmax = $xRange[1]
  $ymin = $yRange[0]
  $ymax = $yRange[1]
  if ($xmin -lt 0) { $xmin = 0.0 }
  if ($ymin -lt 0) { $ymin = 0.0 }

  $axisPen = New-Pen $ink 1.8
  $gridPen = New-Pen $grid 1.0
  $pointBrush = New-Brush $point
  $linePen = New-Pen $line 2.4

  for ($i = 0; $i -le 4; $i++) {
    $gx = $left + (($right - $left) * $i / 4.0)
    $gy = $top + (($bottom - $top) * $i / 4.0)
    $g.DrawLine($gridPen, [single]$gx, [single]$top, [single]$gx, [single]$bottom)
    $g.DrawLine($gridPen, [single]$left, [single]$gy, [single]$right, [single]$gy)
  }

  foreach ($row in $rows) {
    $px = Scale-X ([double]$row.$xMetric) $xmin $xmax $left $right
    $py = Scale-Y ([double]$row.spec_spca_alpha) $ymin $ymax $top $bottom
    $g.FillEllipse($pointBrush, [single]($px - 2.25), [single]($py - 2.25), 4.5, 4.5)
  }

  $stat = $panel.Stat
  $lineX1 = $xmin
  $lineX2 = $xmax
  $lineY1 = [double]$stat.intercept + ([double]$stat.slope * $lineX1)
  $lineY2 = [double]$stat.intercept + ([double]$stat.slope * $lineX2)
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
    Draw-Centered $g (Format-Tick $xVal $xMetric) ($xPos - 42) ($bottom + 6) 84 22 13 $muted
    Draw-Centered $g ("{0:0}" -f $yVal) ($left - 84) ($yPos - 11) 72 22 13 $muted
  }

  $annotation = (Format-R $stat.pearson_r) + "; " + (Format-P $stat.f_p_value)
  Draw-Text $g $annotation ($left + 12) ($top + 10) 15 $ink ([System.Drawing.FontStyle]::Bold)

  $axisPen.Dispose()
  $gridPen.Dispose()
  $pointBrush.Dispose()
  $linePen.Dispose()
}

$dataset = @(Import-Csv -Path (Join-Path $tableDir "spectral_biodiversity_all_metric_dataset.csv"))
$relationships = @(Import-Csv -Path (Join-Path $tableDir "spectral_biodiversity_all_metric_relationships.csv"))

$scales = @("10m", "20m", "50m")
$diversityMetrics = @(
  @{ Name = "phy_faith"; Label = "Faith's PD" },
  @{ Name = "phy_afaith"; Label = "Abundance-weighted Faith's PD" },
  @{ Name = "phy_rao"; Label = "Phylogenetic Rao's Q" },
  @{ Name = "sp_shannon"; Label = "Shannon diversity" },
  @{ Name = "sp_simpson"; Label = "Simpson diversity" }
)

$width = 2800
$height = 3100
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "Vector-normalized alpha-hull area and biodiversity metrics" 150 75 42 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Identity-scale correlations across quadrat grain" 150 135 24 $muted

$lefts = @(430, 1190, 1950)
$tops = @(345, 875, 1405, 1935, 2465)
$panelW = 600
$panelH = 340

for ($c = 0; $c -lt $scales.Count; $c++) {
  Draw-Centered $g ($scales[$c] + " quadrats") $lefts[$c] 235 $panelW 44 28 $ink ([System.Drawing.FontStyle]::Bold)
  $scaleN = (Get-StatRow $relationships $scales[$c] "phy_afaith").n
  Draw-Centered $g ("n = " + $scaleN) $lefts[$c] 275 $panelW 32 20 $muted
}

for ($r = 0; $r -lt $diversityMetrics.Count; $r++) {
  Draw-Text $g $diversityMetrics[$r].Label 430 ($tops[$r] - 36) 23 $ink ([System.Drawing.FontStyle]::Bold)
}

for ($r = 0; $r -lt $diversityMetrics.Count; $r++) {
  for ($c = 0; $c -lt $scales.Count; $c++) {
    $scale = $scales[$c]
    $metric = $diversityMetrics[$r].Name
    $rows = @($dataset | Where-Object {
      $_.scale -eq $scale -and
      $_.spec_spca_alpha -ne "" -and $_.spec_spca_alpha -ne "NA" -and
      $_.$metric -ne "" -and $_.$metric -ne "NA"
    })
    $stat = Get-StatRow $relationships $scale $metric
    Draw-Panel $g @{
      Left = $lefts[$c]
      Top = $tops[$r]
      Right = $lefts[$c] + $panelW
      Bottom = $tops[$r] + $panelH
      Rows = $rows
      DiversityMetric = $metric
      Stat = $stat
    }
  }
}

Draw-Centered $g "Biodiversity metric value" 820 3015 1500 44 27 $ink ([System.Drawing.FontStyle]::Bold)

$font = New-Font 27 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $ink
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("Vector-normalized PCA alpha-hull area", $font, $brush, [System.Drawing.RectangleF]::new(50, 950, 54, 1200), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

$outPath = Join-Path $outDir "11_vector_normalized_alpha_hull_all_diversity_metrics_by_scale.png"
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Output $outPath
