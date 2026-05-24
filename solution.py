from pathlib import Path
from typing import List, Dict, Optional

# --- КОНСТАНТЫ ---
ORIGINAL_FASTQ = Path("reads.fastq")
TRIMMED_SW = Path("trimmed_sw.fastq")
TRIMMED_MINLEN = Path("trimmed_minlen.fastq")


# --- ЗАДАНИЕ 1, 2, 3 ---

def analyze_raw_fastq(filepath: Path) -> None:
    """
    Парсит исходный FASTQ файл за один проход.
    """
    read_lengths = []
    total_gc = 0
    total_bases = 0
    total_phred_pos10 = 0
    read_count_pos10 = 0

    with filepath.open('r', encoding='utf-8') as file:
        for line_num, line in enumerate(file):
            # В FASTQ: строка 1 - ID, 2 - сиквенс (DNA), 3 - "+", 4 - качество (Phred)
            mod = line_num % 4

            # Строка с самой последовательностью ДНК
            if mod == 1:
                sequence = line.strip().upper()
                read_lengths.append(len(sequence))

                # Подсчет для GC-состава (Задание 2)
                total_gc += sequence.count('G') + sequence.count('C')
                total_bases += len(sequence)

            # Строка со значениями качества
            elif mod == 3:
                quality_string = line.strip()
                # Подсчет Phred-качества на 10-й позиции (индекс 9) (Задание 3)
                if len(quality_string) >= 10:
                    # ASCII код символа минус 33 (стандарт Phred+33)
                    phred_score = ord(quality_string[9]) - 33
                    total_phred_pos10 += phred_score
                    read_count_pos10 += 1

    total_reads = len(read_lengths)
    print("--- ЗАДАНИЕ 1 ---")
    print(f"Общее число прочтений: {total_reads}")
    print(f"Минимальная длина: {min(read_lengths)}")
    print(f"Средняя длина: {round(sum(read_lengths) / total_reads)}")
    print(f"Максимальная длина: {max(read_lengths)}\n")

    gc_percentage = (total_gc / total_bases) * 100
    print("--- ЗАДАНИЕ 2 ---")
    print(f"GC-состав: {gc_percentage:.2f}%\n")

    avg_phred = round(total_phred_pos10 / read_count_pos10) if read_count_pos10 else 0
    print("--- ЗАДАНИЕ 3 ---")
    print(f"Среднее качество на 10-й позиции: {avg_phred}\n")


# --- ЗАДАНИЕ 4 ---

def get_read_lengths(filepath: Path) -> List[int]:
    if not filepath.exists():
        return []
    with filepath.open("r", encoding='utf-8') as f:
        # Индекс 1 отвечает за строку с последовательностью
        return [len(line.strip()) for i, line in enumerate(f) if i % 4 == 1]


def get_stats(lengths: List[int]) -> Optional[Dict[str, float]]:
    if not lengths:
        return None
    total_reads = len(lengths)
    return {
        "count": total_reads,
        "min": min(lengths),
        "max": max(lengths),
        "avg_rounded": round(sum(lengths) / total_reads),
    }


def analyze_trimmed_data() -> None:
    orig_lens = get_read_lengths(ORIGINAL_FASTQ)
    sw_lens = get_read_lengths(TRIMMED_SW)
    minlen_lens = get_read_lengths(TRIMMED_MINLEN)

    sw_stats = get_stats(sw_lens)

    print("--- ЗАДАНИЕ 4 ---")
    total_original = len(orig_lens)
    dropped = total_original - sw_stats["count"]

    print(f"Прочтений подверглось триммингу: {dropped}")
    print("Характеристики отфильтрованного файла:")
    print(f" - Минимальная длина: {sw_stats['min']}")
    print(f" - Средняя длина: {sw_stats['avg_rounded']}")
    print(f" - Максимальная длина: {sw_stats['max']}")
    print(f"Тримминг по длине 60, оставшееся число прочтений: {len(minlen_lens)}\n")


def main() -> None:
    analyze_raw_fastq(ORIGINAL_FASTQ)
    analyze_trimmed_data()


if __name__ == "__main__":
    main()
