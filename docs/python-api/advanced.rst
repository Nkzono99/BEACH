実行コンテキストと native kernel
================================================

``Beach.context`` は実行結果と設定の解決を保持します。
``Beach.kernel()`` から受け取る kernel はセッションが所有します。
``field_kernel_from_result`` で作る場合は呼び出し側が閉じます。

.. autoclass:: beach.analysis.context.RunContext
   :members:

.. autofunction:: beach.analysis.field_kernel.field_kernel_from_result

.. autoclass:: beach.native.periodic2.Periodic2Config
   :members:
